//! Optional viral assignments. Derived strain keys never enter public nodes.
use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::Read;
use std::path::Path;
use std::sync::{Arc, Mutex, OnceLock};

use anyhow::{Context, Result};
use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use rayon::prelude::*;
use regex::Regex;
use rusqlite::{Connection, OptionalExtension, params};

use crate::resources;
use crate::{CancellationToken, TaxonId};

const METADATA: &str = "viral.metadata.csv.gz";
const URL: &str =
    "https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/AllNuclMetadata/AllNuclMetadata.csv.gz";
const DB: &str = "nucl.accession2taxid.db";
const BATCH: usize = 16_384;
const MATCH_VERSION: &str = "viral-matching-1";

/// Assignment policy for forward and reverse lookups. Old APIs remain canonical.
#[derive(Clone, Copy, Debug)]
pub struct AccessionMappingOptions {
    pub canonical: bool,
    pub low_memory: bool,
    pub wgs: bool,
    pub threads: Option<usize>,
}
impl Default for AccessionMappingOptions {
    fn default() -> Self {
        Self {
            canonical: true,
            low_memory: true,
            wgs: false,
            threads: None,
        }
    }
}

fn stamp(path: &Path) -> Result<String> {
    let m = fs::metadata(path)?;
    Ok(format!(
        "{}:{}:{:?}",
        path.canonicalize()?.display(),
        m.len(),
        m.modified()?
    ))
}

fn norm(value: &str) -> String {
    value
        .chars()
        .filter(|c| !c.is_whitespace())
        .flat_map(char::to_lowercase)
        .collect()
}
struct Patterns {
    trailing: Regex,
    bare: Regex,
}
impl Patterns {
    fn new() -> Result<Self> {
        Ok(Self {
            trailing: Regex::new(r"\((H\d{1,2}N\d{1,2})[a-z0-9]*\)\)?\s*$")?,
            // Explicit boundary checks reproduce Python lookarounds unsupported by regex.
            bare: Regex::new(r"H[0-9]{1,2}N[0-9]{1,2}")?,
        })
    }
    fn bare<'a>(&self, text: &'a str) -> Option<&'a str> {
        self.bare
            .find_iter(text)
            .find(|m| {
                let before = text[..m.start()].chars().next_back();
                let after = text[m.end()..].chars().next();
                !before.is_some_and(|c| c.is_ascii_alphabetic() || c.is_numeric())
                    && !after.is_some_and(|c| c.is_numeric())
            })
            .map(|m| m.as_str())
    }
    fn subtype<'a>(&self, name: &'a str) -> Option<&'a str> {
        self.trailing
            .captures(name)
            .and_then(|c| c.get(1))
            .map(|m| m.as_str())
            .or_else(|| (!name.contains('(')).then(|| self.bare(name)).flatten())
    }
    fn strain(&self, name: &str) -> Option<String> {
        let inner = name.split_once('(')?.1.trim().strip_suffix(')')?.trim();
        if !inner.starts_with("A/") {
            return None;
        }
        let end = self.trailing.find(inner).map_or(inner.len(), |m| m.start());
        let key = norm(inner[..end].trim());
        (key.len() > 2).then_some(key)
    }
}
struct Matcher {
    strains: HashMap<String, TaxonId>,
    genotypes: HashMap<String, TaxonId>,
    patterns: Patterns,
}
type MatcherCache = Option<(String, Arc<Matcher>)>;
static MATCHER: OnceLock<Mutex<MatcherCache>> = OnceLock::new();
fn matcher(folder: &Path, threads: usize, cancel: &CancellationToken) -> Result<Arc<Matcher>> {
    let key = format!(
        "{}:{}:{}",
        MATCH_VERSION,
        stamp(&folder.join("names.dmp"))?,
        stamp(&folder.join("nodes.dmp"))?
    );
    let cache = MATCHER.get_or_init(|| Mutex::new(None));
    if let Some((old, value)) = cache.lock().unwrap().as_ref() {
        if old == &key {
            return Ok(value.clone());
        }
    }
    cancel.check_cancelled()?;
    let patterns = Patterns::new()?;
    let names = resources::build_names(&folder.join("names.dmp"))?;
    let entries = crate::threads::install(threads, || -> Result<_> {
        let nodes = resources::build_nodes(&folder.join("nodes.dmp"))?;
        nodes
            .par_iter()
            .map(|node| {
                cancel.check_cancelled()?;
                let name = names.get(&node.taxon).map_or("", String::as_str);
                Ok((
                    node.taxon,
                    patterns.strain(name),
                    (node.rank_code == "S4")
                        .then(|| patterns.subtype(name))
                        .flatten()
                        .map(str::to_owned),
                    node.parent,
                ))
            })
            .collect::<Result<Vec<_>>>()
    })??;
    let mut strains = HashMap::new();
    let mut counts: HashMap<String, HashMap<TaxonId, usize>> = HashMap::new();
    for (taxid, strain, subtype, parent) in entries {
        if let Some(strain) = strain {
            strains.entry(strain).or_insert(taxid);
        }
        if let (Some(subtype), Some(parent)) = (subtype, parent) {
            *counts
                .entry(subtype)
                .or_default()
                .entry(parent)
                .or_default() += 1;
        }
    }
    let genotypes = counts
        .into_iter()
        .map(|(subtype, parents)| {
            let parent = parents
                .into_iter()
                .max_by(|(a, n), (b, m)| n.cmp(m).then_with(|| b.cmp(a)))
                .unwrap()
                .0;
            (subtype, parent)
        })
        .collect();
    let value = Arc::new(Matcher {
        strains,
        genotypes,
        patterns,
    });
    *cache.lock().unwrap() = Some((key, value.clone()));
    Ok(value)
}

fn columns(reader: &mut csv::Reader<impl Read>) -> Result<[usize; 3]> {
    let header = reader.headers()?;
    let mut columns = [0; 3];
    for (slot, name) in columns.iter_mut().zip(["#Accession", "Genotype", "Strain"]) {
        *slot = header
            .iter()
            .position(|s| s == name)
            .with_context(|| format!("viral metadata missing {name}"))?;
    }
    Ok(columns)
}
fn present(text: &str) -> bool {
    !text.trim().is_empty()
}
fn trim_metadata(
    input: impl Read,
    output: impl std::io::Write,
    cancel: &CancellationToken,
) -> Result<()> {
    let mut reader = csv::Reader::from_reader(input);
    let [a, g, s] = columns(&mut reader)?;
    let mut writer = csv::Writer::from_writer(output);
    writer.write_record(["#Accession", "Genotype", "Strain"])?;
    for (i, row) in reader.records().enumerate() {
        if i % BATCH == 0 {
            cancel.check_cancelled()?;
        }
        let row = row?;
        if present(&row[a]) && (present(&row[g]) || present(&row[s])) {
            writer.write_record([&row[a], &row[g], &row[s]])?;
        }
    }
    writer.flush()?;
    cancel.check_cancelled()
}
fn metadata(folder: &Path, refresh: bool, cancel: &CancellationToken) -> Result<()> {
    let path = folder.join(METADATA);
    if path.exists() && !refresh {
        return Ok(());
    }
    fs::create_dir_all(folder)?;
    cancel.check_cancelled()?;
    let response = reqwest::blocking::Client::builder()
        .timeout(std::time::Duration::from_secs(1800))
        .build()?
        .get(URL)
        .send()?
        .error_for_status()?;
    install_metadata(folder, MultiGzDecoder::new(response), cancel)
}

fn install_metadata(folder: &Path, input: impl Read, cancel: &CancellationToken) -> Result<()> {
    cancel.check_cancelled()?;
    let mut temp = tempfile::NamedTempFile::new_in(folder)?;
    let mut gzip = GzEncoder::new(temp.as_file_mut(), Compression::default());
    trim_metadata(input, &mut gzip, cancel)?;
    gzip.finish()?;
    cancel.check_cancelled()?;
    temp.persist(folder.join(METADATA))?;
    Ok(())
}

// Priority 0 is a strain match; priority 1 is genotype. Ordered reduction makes
// duplicate handling independent of worker scheduling and batch boundaries.
fn scan(
    folder: &Path,
    requested: Option<&HashSet<String>>,
    threads: usize,
    cancel: &CancellationToken,
    mut visit: impl FnMut(&str, TaxonId, u8) -> Result<()>,
) -> Result<()> {
    let matcher = matcher(folder, threads, cancel)?;
    let decoder = resources::lookup_decoder(threads)?;
    let mut reader = csv::Reader::from_reader(decoder.open(folder.join(METADATA))?);
    let [a, g, s] = columns(&mut reader)?;
    let mut rows = reader.records();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?;
    loop {
        cancel.check_cancelled()?;
        let batch = rows
            .by_ref()
            .take(BATCH)
            .collect::<std::result::Result<Vec<_>, _>>()?;
        if batch.is_empty() {
            break;
        }
        let hits = pool.install(|| {
            batch
                .par_iter()
                .map(|row| {
                    if requested.is_some_and(|wanted| !wanted.contains(&row[a])) {
                        return None;
                    }
                    matcher
                        .strains
                        .get(&norm(&row[s]))
                        .map(|t| (*t, 0))
                        .or_else(|| {
                            matcher
                                .patterns
                                .bare(&row[g])
                                .and_then(|g| matcher.genotypes.get(g))
                                .map(|t| (*t, 1))
                        })
                })
                .collect::<Vec<_>>()
        });
        for (row, hit) in batch.iter().zip(hits) {
            if let Some((taxid, priority)) = hit {
                visit(&row[a], taxid, priority)?;
            }
        }
    }
    cancel.check_cancelled()
}
fn winners(
    folder: &Path,
    requested: Option<&HashSet<String>>,
    threads: usize,
    cancel: &CancellationToken,
) -> Result<HashMap<String, (TaxonId, u8)>> {
    let mut found: HashMap<String, (TaxonId, u8)> = HashMap::new();
    scan(folder, requested, threads, cancel, |a, t, p| {
        let entry = found.entry(a.to_owned()).or_insert((t, p));
        if p < entry.1 {
            *entry = (t, p);
        }
        Ok(())
    })?;
    Ok(found)
}

fn fingerprint(folder: &Path, connection: &Connection) -> Result<String> {
    let mut key = MATCH_VERSION.to_owned();
    for file in ["names.dmp", "nodes.dmp", METADATA] {
        key.push_str(&stamp(&folder.join(file))?);
    }
    for file in ["nucl_gb.accession2taxid.gz", "nucl_wgs.accession2taxid.gz"] {
        if folder.join(file).exists() {
            key.push_str(&stamp(&folder.join(file))?);
        }
    }
    let mut stmt = connection.prepare("SELECT source, status, etag, last_modified, size, row_count FROM a2t_sources ORDER BY source")?;
    let mut rows = stmt.query([])?;
    while let Some(row) = rows.next()? {
        for i in 0..6 {
            key.push_str(&format!("{:?};", row.get_ref(i)?));
        }
    }
    Ok(key)
}
fn ensure_overrides(
    folder: &Path,
    wgs: bool,
    threads: Option<usize>,
    cancel: &CancellationToken,
) -> Result<()> {
    resources::ensure_accession_database_with_cancel(
        folder,
        crate::AccessionDatabaseOptions {
            wgs,
            threads,
            ..Default::default()
        },
        cancel,
    )?;
    let mut connection = Connection::open(folder.join(DB))?;
    connection.busy_timeout(std::time::Duration::from_secs(60))?;
    let token = cancel.clone();
    connection.progress_handler(10_000, Some(move || token.is_cancelled()))?;
    // Serialize checking/building so another caller never observes a partial generation.
    let transaction =
        connection.transaction_with_behavior(rusqlite::TransactionBehavior::Immediate)?;
    let key = fingerprint(folder, &transaction)?;
    let old: Option<String> = transaction
        .query_row(
            "SELECT value FROM a2t_meta WHERE key='viral_overrides'",
            [],
            |r| r.get(0),
        )
        .optional()?;
    if old.as_deref() == Some(&key) {
        return Ok(());
    }
    transaction.execute_batch("CREATE TABLE IF NOT EXISTS a2t_overrides (accession TEXT PRIMARY KEY, taxid INTEGER NOT NULL) WITHOUT ROWID;
        CREATE INDEX IF NOT EXISTS idx_override_taxid ON a2t_overrides(taxid);
        CREATE TEMP TABLE viral_candidates (accession TEXT PRIMARY KEY, taxid INTEGER NOT NULL, priority INTEGER NOT NULL) WITHOUT ROWID;")?;
    {
        let mut insert = transaction.prepare("INSERT INTO viral_candidates VALUES (?1,?2,?3) ON CONFLICT(accession) DO UPDATE SET taxid=excluded.taxid, priority=excluded.priority WHERE excluded.priority < viral_candidates.priority")?;
        scan(
            folder,
            None,
            crate::threads::resolve(threads)?,
            cancel,
            |a, t, p| {
                insert.execute(params![a, t, p])?;
                Ok(())
            },
        )?;
    }
    cancel.check_cancelled()?;
    transaction.execute_batch("DELETE FROM a2t_overrides;
        INSERT INTO a2t_overrides SELECT v.accession,v.taxid FROM viral_candidates v LEFT JOIN a2t a ON a.accession=v.accession WHERE a.taxid IS NULL OR a.taxid != v.taxid;
        DROP TABLE viral_candidates;")?;
    transaction.execute(
        "INSERT OR REPLACE INTO a2t_meta VALUES ('viral_overrides', ?1)",
        [&key],
    )?;
    cancel.check_cancelled()?;
    transaction.commit()?;
    Ok(())
}

/// Prepare optional resources. Canonical users only refresh an already-installed
/// metadata cache; they never download it for the first time.
pub fn prepare_alternative_mappings(
    folder: impl AsRef<Path>,
    canonical: bool,
    low_memory: bool,
    wgs: bool,
    refresh: bool,
    threads: Option<usize>,
    cancel: &CancellationToken,
) -> Result<()> {
    let folder = folder.as_ref();
    if canonical {
        if refresh && folder.join(METADATA).exists() {
            metadata(folder, true, cancel)?;
        }
        return Ok(());
    }
    metadata(folder, refresh, cancel).context("cannot prepare alternative viral mappings")?;
    if !low_memory {
        ensure_overrides(folder, wgs, threads, cancel)?;
    }
    Ok(())
}

pub fn lookup_accession_taxids_with_options(
    folder: impl AsRef<Path>,
    requested: HashSet<String>,
    options: AccessionMappingOptions,
    cancel: &CancellationToken,
) -> Result<HashMap<String, TaxonId>> {
    if requested.is_empty() {
        return Ok(HashMap::new());
    }
    let folder = folder.as_ref();
    let mut found = resources::lookup_accession_taxids_with_cancel(
        folder,
        requested.clone(),
        options.low_memory,
        options.wgs,
        options.threads,
        cancel,
    )?;
    if options.canonical {
        return Ok(found);
    }
    prepare_alternative_mappings(
        folder,
        false,
        options.low_memory,
        options.wgs,
        false,
        options.threads,
        cancel,
    )?;
    if options.low_memory {
        for (a, (t, _)) in winners(
            folder,
            Some(&requested),
            crate::threads::resolve(options.threads)?,
            cancel,
        )? {
            found.insert(a, t);
        }
    } else {
        let mut connection = Connection::open(folder.join(DB))?;
        let tx = connection.transaction()?;
        tx.execute_batch("CREATE TEMP TABLE wanted (accession TEXT PRIMARY KEY) WITHOUT ROWID;")?;
        {
            let mut insert = tx.prepare("INSERT INTO wanted VALUES (?1)")?;
            for a in requested {
                insert.execute([a])?;
            }
        }
        let mut stmt = tx.prepare("SELECT o.accession,o.taxid FROM wanted w CROSS JOIN a2t_overrides o ON o.accession=w.accession")?;
        let mut rows = stmt.query([])?;
        while let Some(row) = rows.next()? {
            cancel.check_cancelled()?;
            found.insert(row.get(0)?, row.get(1)?);
        }
    }
    Ok(found)
}

pub fn lookup_taxid_accessions_with_options(
    folder: impl AsRef<Path>,
    taxa: &[TaxonId],
    options: AccessionMappingOptions,
    cancel: &CancellationToken,
) -> Result<HashSet<String>> {
    if taxa.is_empty() {
        return Ok(HashSet::new());
    }
    let folder = folder.as_ref();
    if options.canonical {
        return resources::lookup_taxid_accessions_with_cancel(
            folder,
            taxa,
            options.low_memory,
            options.wgs,
            options.threads,
            cancel,
        );
    }
    prepare_alternative_mappings(
        folder,
        false,
        options.low_memory,
        options.wgs,
        false,
        options.threads,
        cancel,
    )?;
    if options.low_memory {
        let mut found = resources::lookup_taxid_accessions_with_cancel(
            folder,
            taxa,
            true,
            options.wgs,
            options.threads,
            cancel,
        )?;
        let wanted: HashSet<_> = taxa.iter().copied().collect();
        let threads = crate::threads::resolve(options.threads)?;
        // First discover possible results, then resolve all duplicate priorities
        // for just those accessions. Two scans avoid retaining every viral override
        // when a reverse query requests only a small taxonomic subset.
        scan(folder, None, threads, cancel, |a, t, _| {
            if wanted.contains(&t) {
                found.insert(a.to_owned());
            }
            Ok(())
        })?;
        if !found.is_empty() {
            for (a, (t, _)) in winners(folder, Some(&found), threads, cancel)? {
                if !wanted.contains(&t) {
                    found.remove(&a);
                }
            }
        }
        return Ok(found);
    }
    let mut connection = Connection::open(folder.join(DB))?;
    let tx = connection.transaction()?;
    tx.execute_batch("CREATE TEMP TABLE wanted_taxa (taxid INTEGER PRIMARY KEY);")?;
    {
        let mut insert = tx.prepare("INSERT OR IGNORE INTO wanted_taxa VALUES (?1)")?;
        for t in taxa {
            insert.execute([t])?;
        }
    }
    let mut stmt = tx.prepare("SELECT a.accession FROM wanted_taxa w CROSS JOIN a2t a INDEXED BY idx_taxid ON a.taxid=w.taxid
        WHERE NOT EXISTS (SELECT 1 FROM a2t_overrides o WHERE o.accession=a.accession)
        UNION ALL SELECT o.accession FROM wanted_taxa w CROSS JOIN a2t_overrides o INDEXED BY idx_override_taxid ON o.taxid=w.taxid")?;
    let mut rows = stmt.query([])?;
    let mut found = HashSet::new();
    while let Some(row) = rows.next()? {
        cancel.check_cancelled()?;
        found.insert(row.get(0)?);
    }
    Ok(found)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn gzip(path: &Path, text: &str) {
        let mut out = GzEncoder::new(std::fs::File::create(path).unwrap(), Compression::fast());
        out.write_all(text.as_bytes()).unwrap();
        out.finish().unwrap();
    }
    fn fixture() -> tempfile::TempDir {
        let dir = tempfile::tempdir().unwrap();
        fs::write(dir.path().join("nodes.dmp"), "1 | 1 | no rank |\n10 | 1 | species |\n11 | 10 | no rank |\n12 | 11 | no rank |\n13 | 12 | no rank |\n14 | 11 | no rank |\n15 | 14 | no rank |\n16 | 12 | no rank |\n").unwrap();
        fs::write(dir.path().join("names.dmp"), "1 | root | | scientific name |\n10 | Influenza A virus | | scientific name |\n11 | intermediate | | scientific name |\n12 | subtype parent | | scientific name |\n13 | Influenza A virus (A/Test/1/2020(H1N1)) | | scientific name |\n14 | other parent | | scientific name |\n15 | Influenza A virus (A/Other/2/2020(H1N1)) | | scientific name |\n16 | Influenza A virus (A/Test/1/2020(H2N2)) | | scientific name |\n").unwrap();
        fs::write(dir.path().join("targets.json"), "{\"pathogens\":{}}").unwrap();
        gzip(
            &dir.path().join("nucl_gb.accession2taxid.gz"),
            "accession\taccession.version\ttaxid\tgi\nNC_000001\tNC_000001.1\t10\t0\nNC_000002\tNC_000002.1\t10\t0\nNC_000003\tNC_000003.1\t10\t0\nNC_000004\tNC_000004.1\t10\t0\nNC_000005\tNC_000005.1\t10\t0\nNC_000007\tNC_000007.1\t13\t0\n",
        );
        gzip(
            &dir.path().join(METADATA),
            "#Accession,Genotype,Strain\nNC_000001.1,H1N1, A/test/1/2020 \nNC_000002.1,H1N1,unknown\nNC_000003.1,,A/Test/1/2020\nNC_000004.1,H99N99,unknown\nNC_000005.1,H1N1,unknown\nNC_000005.1,H1N1,A/Test/1/2020\nNC_000005.1,H1N1,A/Other/2/2020\nNC_000006.1,,A/Test/1/2020\nNC_000007.1,,A/Test/1/2020\n",
        );
        dir
    }
    fn requests() -> HashSet<String> {
        (1..=7).map(|i| format!("NC_{i:06}.1")).collect()
    }
    #[test]
    fn parsing_boundaries_and_private_strains() {
        let p = Patterns::new().unwrap();
        assert_eq!(
            p.strain("Influenza A virus (A/Test/1/2020(H1N1))")
                .as_deref(),
            Some("a/test/1/2020")
        );
        assert_eq!(p.strain("Influenza A virus (A/)"), None);
        assert_eq!(p.strain("Influenza B virus (B/Test/1/2020)"), None);
        assert_eq!(
            p.subtype("Influenza A virus (A/Test(H1N1pdm09))"),
            Some("H1N1")
        );
        assert_eq!(p.subtype("Influenza A H1N1"), Some("H1N1"));
        assert_eq!(p.subtype("Influenza A (H1N1 extra)"), None);
        assert_eq!(p.bare("AH1N1"), None);
        assert_eq!(p.bare("H1N123"), None);
        assert_eq!(p.bare("H1N1pdm09"), Some("H1N1"));
    }
    #[test]
    fn modes_threads_priorities_and_reverse_agree() {
        let dir = fixture();
        let expected: HashMap<_, _> = [13, 12, 13, 10, 13, 13, 13]
            .into_iter()
            .enumerate()
            .map(|(i, t)| (format!("NC_{:06}.1", i + 1), t))
            .collect();
        for low_memory in [true, false] {
            for threads in [1, 3] {
                let options = AccessionMappingOptions {
                    canonical: false,
                    low_memory,
                    threads: Some(threads),
                    ..Default::default()
                };
                let found = lookup_accession_taxids_with_options(
                    dir.path(),
                    requests(),
                    options,
                    &CancellationToken::default(),
                )
                .unwrap();
                assert_eq!(found, expected);
                for t in [10, 12, 13, 14, 15, 999] {
                    assert_eq!(
                        lookup_taxid_accessions_with_options(
                            dir.path(),
                            &[t],
                            options,
                            &CancellationToken::default()
                        )
                        .unwrap(),
                        expected
                            .iter()
                            .filter(|(_, v)| **v == t)
                            .map(|(a, _)| a.clone())
                            .collect()
                    );
                }
            }
            if low_memory {
                assert!(!dir.path().join(DB).exists());
            }
        }
        let c = Connection::open(dir.path().join(DB)).unwrap();
        assert_eq!(
            c.query_row("SELECT count(*) FROM a2t_overrides", [], |r| r
                .get::<_, i64>(0))
                .unwrap(),
            5
        );
        assert_eq!(
            c.query_row(
                "SELECT count(*) FROM a2t_overrides WHERE accession='NC_000007.1'",
                [],
                |r| r.get::<_, i64>(0)
            )
            .unwrap(),
            0
        );
    }
    #[test]
    fn canonical_needs_no_metadata_and_rust_object_uses_policy() {
        let dir = fixture();
        fs::remove_file(dir.path().join(METADATA)).unwrap();
        let found = lookup_accession_taxids_with_options(
            dir.path(),
            requests(),
            Default::default(),
            &CancellationToken::default(),
        )
        .unwrap();
        assert_eq!(found["NC_000001.1"], 10);
        assert!(!dir.path().join(METADATA).exists());
        let dir = fixture();
        let mut tu = crate::TaxutilsBuilder::new()
            .save_folder(dir.path())
            .canonical(false)
            .threads(Some(2))
            .accessions(["NC_000001.1".to_owned()])
            .build()
            .unwrap();
        assert_eq!(tu.a2t["NC_000001.1"], 13);
        tu.load_a2t(&["NC_000002.1"], None, true, None).unwrap();
        assert_eq!(tu.a2t.len(), 2);
        assert_eq!(tu.a2t["NC_000002.1"], 12);
        assert!(
            tu.get_t2a(&[13], None, None)
                .unwrap()
                .contains("NC_000001.1")
        );
    }
    #[test]
    fn replacing_metadata_removes_obsolete_overrides_and_failure_rolls_back() {
        let dir = fixture();
        let options = AccessionMappingOptions {
            canonical: false,
            low_memory: false,
            threads: Some(2),
            ..Default::default()
        };
        lookup_accession_taxids_with_options(
            dir.path(),
            requests(),
            options,
            &CancellationToken::default(),
        )
        .unwrap();
        gzip(
            &dir.path().join(METADATA),
            "#Accession,Genotype,Strain\nNC_000002.1,H1N1,unknown\n",
        );
        let found = lookup_accession_taxids_with_options(
            dir.path(),
            requests(),
            options,
            &CancellationToken::default(),
        )
        .unwrap();
        assert_eq!(found["NC_000001.1"], 10);
        assert_eq!(found["NC_000002.1"], 12);
        assert!(!found.contains_key("NC_000006.1"));
        gzip(&dir.path().join(METADATA), "wrong,columns\na,b\n");
        assert!(
            lookup_accession_taxids_with_options(
                dir.path(),
                requests(),
                options,
                &CancellationToken::default()
            )
            .is_err()
        );
        let c = Connection::open(dir.path().join(DB)).unwrap();
        assert_eq!(
            c.query_row("SELECT count(*) FROM a2t_overrides", [], |r| r
                .get::<_, i64>(0))
                .unwrap(),
            1
        );
    }
    #[test]
    fn trimming_handles_csv_quoting_and_missing_genotype() {
        let mut output = Vec::new();
        trim_metadata(
            "#Accession,Genotype,Strain,Unused\nNC_1.1,,\"A/Test, place/1\",x\nNC_2.1,,,y\n"
                .as_bytes(),
            &mut output,
            &CancellationToken::default(),
        )
        .unwrap();
        let rows: Vec<_> = csv::Reader::from_reader(output.as_slice())
            .records()
            .map(Result::unwrap)
            .collect();
        assert_eq!(rows.len(), 1);
        assert_eq!(&rows[0][2], "A/Test, place/1");
        let cancel = CancellationToken::default();
        cancel.cancel();
        assert!(
            trim_metadata(
                "#Accession,Genotype,Strain\n".as_bytes(),
                Vec::new(),
                &cancel
            )
            .is_err()
        );
    }
    #[test]
    fn failed_metadata_install_and_cancel_preserve_cache() {
        let dir = fixture();
        let original = fs::read(dir.path().join(METADATA)).unwrap();
        assert!(
            install_metadata(
                dir.path(),
                "bad,header\na,b\n".as_bytes(),
                &CancellationToken::default()
            )
            .is_err()
        );
        assert_eq!(fs::read(dir.path().join(METADATA)).unwrap(), original);
        let cancel = CancellationToken::default();
        cancel.cancel();
        assert!(
            install_metadata(
                dir.path(),
                "#Accession,Genotype,Strain\n".as_bytes(),
                &cancel
            )
            .is_err()
        );
        assert_eq!(fs::read(dir.path().join(METADATA)).unwrap(), original);
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 5);
    }

    #[test]
    fn taxonomy_and_canonical_source_changes_invalidate_overrides() {
        let dir = fixture();
        let options = AccessionMappingOptions {
            canonical: false,
            low_memory: false,
            threads: Some(1),
            ..Default::default()
        };
        let cancel = CancellationToken::default();
        lookup_accession_taxids_with_options(dir.path(), requests(), options, &cancel).unwrap();
        // Simulate a committed canonical refresh, including its source validator.
        let connection = Connection::open(dir.path().join(DB)).unwrap();
        connection.execute_batch("BEGIN; UPDATE a2t SET taxid=10 WHERE accession='NC_000007.1'; UPDATE a2t_sources SET etag='new-generation'; COMMIT;").unwrap();
        let found =
            lookup_accession_taxids_with_options(dir.path(), requests(), options, &cancel).unwrap();
        assert_eq!(found["NC_000007.1"], 13); // previously redundant, now an override
        let path = dir.path().join("names.dmp");
        fs::write(
            &path,
            fs::read_to_string(&path)
                .unwrap()
                .replace("A/Test/1/2020", "A/Changed/1/2020"),
        )
        .unwrap();
        let found =
            lookup_accession_taxids_with_options(dir.path(), requests(), options, &cancel).unwrap();
        assert_eq!(found["NC_000003.1"], 10); // missing genotype, old strain is gone
        assert_eq!(found["NC_000001.1"], 12); // now genotype fallback
        assert!(!found.contains_key("NC_000006.1"));
    }

    #[test]
    fn duplicate_priority_crosses_batch_boundary_and_wgs_upgrade() {
        let dir = fixture();
        let mut metadata = String::from("#Accession,Genotype,Strain\nNC_000001.1,H1N1,unknown\n");
        for _ in 0..BATCH {
            metadata.push_str("NC_999999.1,,unknown\n");
        }
        metadata.push_str("NC_000001.1,,A/Test/1/2020\nABCD01000001.1,,A/Test/1/2020\n");
        gzip(&dir.path().join(METADATA), &metadata);
        gzip(
            &dir.path().join("nucl_wgs.accession2taxid.gz"),
            "accession\taccession.version\ttaxid\tgi\nABCD01000001\tABCD01000001.1\t10\t0\n",
        );
        let mut requested = requests();
        requested.insert("ABCD01000001.1".to_owned());
        for low_memory in [true, false] {
            let options = AccessionMappingOptions {
                canonical: false,
                low_memory,
                wgs: true,
                threads: Some(3),
            };
            let found = lookup_accession_taxids_with_options(
                dir.path(),
                requested.clone(),
                options,
                &CancellationToken::default(),
            )
            .unwrap();
            assert_eq!(found["NC_000001.1"], 13);
            assert_eq!(found["ABCD01000001.1"], 13);
            assert!(
                !lookup_taxid_accessions_with_options(
                    dir.path(),
                    &[12],
                    options,
                    &CancellationToken::default()
                )
                .unwrap()
                .contains("NC_000001.1")
            );
        }
    }
}
