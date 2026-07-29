use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::read::GzDecoder;
use rayon::prelude::*;
use reqwest::blocking::Client;
use rusqlite::{Connection, params};
use serde_json::Value;

use crate::accession::parse_accessions;
use crate::taxonomy::{
    TaxonId, TaxonNode, TaxonomicUtils, assign_rank_codes, canonical_name, rank_index,
    taxonomic_order,
};

const TAXDUMP_URL: &str = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";
const TARGETS_URL: &str = "https://web.cs.ucla.edu/~wob/projects/taxutils/targets.json";
const A2T_BASE_URL: &str = "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid";
const GB_FILE: &str = "nucl_gb.accession2taxid.gz";
const WGS_FILE: &str = "nucl_wgs.accession2taxid.gz";
const DB_FILE: &str = "nucl.accession2taxid.db";

#[derive(Clone, Debug)]
pub struct TaxutilsOptions {
    pub accessions: Option<Vec<String>>,
    pub low_memory: bool,
    pub targets_json: Option<PathBuf>,
    pub rebuild: bool,
    pub wgs: bool,
    pub save_folder: PathBuf,
}

impl Default for TaxutilsOptions {
    fn default() -> Self {
        Self {
            accessions: None,
            low_memory: true,
            targets_json: None,
            rebuild: false,
            wgs: false,
            save_folder: std::env::var_os("TAXUTILS_GLOBALS")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from("./taxutils/")),
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct TaxutilsBuilder {
    options: TaxutilsOptions,
}

impl TaxutilsBuilder {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn accessions(mut self, accessions: impl IntoIterator<Item = String>) -> Self {
        self.options.accessions = Some(accessions.into_iter().collect());
        self
    }
    pub fn low_memory(mut self, value: bool) -> Self {
        self.options.low_memory = value;
        self
    }
    pub fn targets_json(mut self, value: impl Into<PathBuf>) -> Self {
        self.options.targets_json = Some(value.into());
        self
    }
    pub fn rebuild(mut self, value: bool) -> Self {
        self.options.rebuild = value;
        self
    }
    pub fn wgs(mut self, value: bool) -> Self {
        self.options.wgs = value;
        self
    }
    pub fn save_folder(mut self, value: impl Into<PathBuf>) -> Self {
        self.options.save_folder = value.into();
        self
    }
    pub fn build(self) -> Result<TaxonomicUtils> {
        TaxonomicUtils::new(self.options)
    }
}

pub(crate) fn load_taxutils(options: TaxutilsOptions) -> Result<TaxonomicUtils> {
    fs::create_dir_all(&options.save_folder)?;
    let names_path = options.save_folder.join("names.dmp");
    let nodes_path = options.save_folder.join("nodes.dmp");
    if options.rebuild || !names_path.exists() || !nodes_path.exists() {
        download_taxonomy(&options.save_folder, &names_path, &nodes_path)?;
    }
    let targets_path = if let Some(path) = &options.targets_json {
        path.clone()
    } else {
        let path = options.save_folder.join("targets.json");
        if options.rebuild || !path.exists() {
            download_file(TARGETS_URL, &path)?;
        }
        path
    };
    let mut names = build_names(&names_path)?;
    let nodes = build_nodes(&nodes_path)?;
    let target_taxa = build_target_taxa(&nodes, &names, &targets_path)?;
    if options.rebuild || !options.low_memory {
        ensure_a2t_db(&options.save_folder, options.rebuild, options.wgs)?;
    }
    let mut a2t = HashMap::new();
    if let Some(accessions) = &options.accessions {
        a2t = lookup_a2t(
            &options.save_folder,
            accessions,
            options.low_memory,
            options.wgs,
        )?;
    }
    names.insert(2697049, "SARS-CoV-2".to_owned());
    names.insert(694009, "SARS-related-CoV".to_owned());
    Ok(TaxonomicUtils::from_parts(
        names,
        nodes,
        target_taxa,
        a2t,
        options.low_memory,
        options.wgs,
        options.save_folder,
    ))
}

impl TaxonomicUtils {
    /// Load accession-to-taxid mappings. By default this replaces `a2t`.
    pub fn load_a2t<S: AsRef<str> + Sync>(
        &mut self,
        accessions: &[S],
        low_memory: Option<bool>,
        extend: bool,
        wgs: Option<bool>,
    ) -> Result<()> {
        let mut requested = parse_accessions(accessions, true)
            .into_iter()
            .filter(|value| value != "NA")
            .collect::<HashSet<_>>();
        if extend {
            requested.retain(|accession| !self.a2t.contains_key(accession));
            if requested.is_empty() {
                return Ok(());
            }
        }
        let found = lookup_a2t(
            &self.save_folder,
            &requested.into_iter().collect::<Vec<_>>(),
            low_memory.unwrap_or(self.low_memory),
            wgs.unwrap_or(self.wgs),
        )?;
        if !extend {
            self.a2t.clear();
        }
        self.a2t.extend(found);
        Ok(())
    }

    pub fn get_t2a(
        &self,
        taxa: &[TaxonId],
        low_memory: Option<bool>,
        wgs: Option<bool>,
    ) -> Result<HashSet<String>> {
        lookup_t2a(
            &self.save_folder,
            taxa,
            low_memory.unwrap_or(self.low_memory),
            wgs.unwrap_or(self.wgs),
        )
    }
}

fn download_file(url: &str, path: &Path) -> Result<()> {
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let mut temporary = tempfile::NamedTempFile::new_in(parent)?;
    let mut response = Client::new()
        .get(url)
        .send()
        .with_context(|| format!("failed to download {url}"))?
        .error_for_status()?;
    std::io::copy(&mut response, &mut temporary)?;
    temporary.flush()?;
    temporary
        .persist(path)
        .map_err(|error| error.error)
        .with_context(|| format!("failed to install {}", path.display()))?;
    Ok(())
}

fn download_taxonomy(save_folder: &Path, names_path: &Path, nodes_path: &Path) -> Result<()> {
    let tarball = save_folder.join("taxdump.tar.gz");
    download_file(TAXDUMP_URL, &tarball)?;
    let file = File::open(&tarball)?;
    let mut archive = tar::Archive::new(GzDecoder::new(file));
    let mut found = HashSet::new();
    for entry in archive.entries()? {
        let mut entry = entry?;
        let path = entry.path()?;
        let Some(filename) = path
            .file_name()
            .and_then(|value| value.to_str())
            .map(str::to_owned)
        else {
            continue;
        };
        let destination = match filename.as_str() {
            "names.dmp" => names_path,
            "nodes.dmp" => nodes_path,
            _ => continue,
        };
        let mut output = File::create(destination)?;
        std::io::copy(&mut entry, &mut output)?;
        found.insert(filename);
    }
    let _ = fs::remove_file(&tarball);
    if !found.contains("names.dmp") || !found.contains("nodes.dmp") {
        bail!("could not find names.dmp and nodes.dmp in taxdump")
    }
    Ok(())
}

fn build_names(path: &Path) -> Result<HashMap<TaxonId, String>> {
    let mut names = HashMap::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let mut fields = line.split('|').map(str::trim);
        let (Some(taxon), Some(name), Some(_unique_name), Some(name_class)) =
            (fields.next(), fields.next(), fields.next(), fields.next())
        else {
            continue;
        };
        if name_class == "scientific name" {
            names.insert(taxon.parse()?, name.to_owned());
        }
    }
    names.insert(0, "unclassified".to_owned());
    Ok(names)
}

fn build_nodes(path: &Path) -> Result<Vec<TaxonNode>> {
    let mut raw = Vec::new();
    let mut parent = HashMap::new();
    let mut ranks = HashMap::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let mut fields = line.split('|').map(str::trim);
        let (Some(taxon), Some(raw_parent), Some(rank)) =
            (fields.next(), fields.next(), fields.next())
        else {
            continue;
        };
        let taxon: TaxonId = taxon.parse()?;
        let raw_parent: TaxonId = raw_parent.parse()?;
        let parent_taxon = Some(raw_parent);
        let rank = rank.to_ascii_lowercase();
        raw.push((taxon, parent_taxon, rank.clone()));
        parent.insert(taxon, parent_taxon);
        ranks.insert(taxon, rank);
    }
    let codes = assign_rank_codes(&parent, &ranks);
    Ok(raw
        .into_par_iter()
        .map(|(taxon, parent, rank)| {
            let rank_code = codes[&taxon].clone();
            let rank_base = rank_code.chars().next().unwrap_or('U');
            TaxonNode {
                taxon,
                parent,
                rank,
                rank_code,
                rank_base,
                rank_idx: rank_index(rank_base),
                new_rank: canonical_name(rank_base).to_owned(),
            }
        })
        .collect())
}

fn build_target_taxa(
    nodes: &[TaxonNode],
    names: &HashMap<TaxonId, String>,
    path: &Path,
) -> Result<Vec<TaxonId>> {
    let value: Value = serde_json::from_reader(File::open(path)?)?;
    let pathogens = value
        .get("pathogens")
        .and_then(Value::as_object)
        .context("targets JSON has no pathogens object")?;
    let pathogen_taxa = pathogens
        .values()
        .map(|value| {
            value
                .as_i64()
                .or_else(|| value.as_str().and_then(|text| text.parse().ok()))
                .context("invalid pathogen taxid")
        })
        .collect::<Result<Vec<_>>>()?;
    let mut parent = nodes
        .iter()
        .map(|node| (node.taxon, node.parent))
        .collect::<HashMap<_, _>>();
    parent.insert(1, None);
    let mut children: HashMap<TaxonId, Vec<TaxonId>> = HashMap::new();
    for node in nodes {
        if node.taxon != 1
            && let Some(parent) = parent.get(&node.taxon).copied().flatten()
        {
            children.entry(parent).or_default().push(node.taxon);
        }
    }
    let rank_idx = nodes
        .iter()
        .map(|node| (node.taxon, node.rank_idx))
        .collect::<HashMap<_, _>>();
    let ranks = nodes
        .iter()
        .map(|node| (node.taxon, node.rank_code.clone()))
        .collect();
    let taxa = pathogen_taxa
        .par_iter()
        .map(|pathogen| {
            let mut local = HashSet::new();
            let mut stack = vec![*pathogen];
            while let Some(node) = stack.pop() {
                local.insert(node);
                if let Some(values) = children.get(&node) {
                    stack.extend(values.iter().rev().copied());
                }
            }
            let mut current = parent.get(pathogen).copied().flatten();
            while let Some(node) = current {
                if rank_idx.get(&node).copied().unwrap_or(6) < 7 {
                    break;
                }
                local.insert(node);
                current = parent.get(&node).copied().flatten();
            }
            local
        })
        .reduce(HashSet::new, |mut all, local| {
            all.extend(local);
            all
        });
    Ok(taxonomic_order(taxa, &parent, &ranks, names))
}

fn a2t_paths(save_folder: &Path, wgs: bool) -> Vec<PathBuf> {
    let mut paths = vec![save_folder.join(GB_FILE)];
    if wgs {
        paths.push(save_folder.join(WGS_FILE));
    }
    paths
}
fn ensure_a2t_files(save_folder: &Path, wgs: bool, rebuild: bool) -> Result<Vec<PathBuf>> {
    let paths = a2t_paths(save_folder, wgs);
    paths.par_iter().try_for_each(|path| {
        if rebuild || !path.exists() {
            let filename = path.file_name().unwrap().to_string_lossy();
            download_file(&format!("{A2T_BASE_URL}/{filename}"), path)?;
        }
        Ok::<_, anyhow::Error>(())
    })?;
    Ok(paths)
}
fn a2t_columns<R: BufRead>(reader: &mut R) -> Result<(usize, usize)> {
    let mut header = String::new();
    reader.read_line(&mut header)?;
    let accession = header
        .trim_end()
        .split('\t')
        .position(|value| value == "accession.version")
        .context("missing accession.version column")?;
    let taxid = header
        .trim_end()
        .split('\t')
        .position(|value| value == "taxid")
        .context("missing taxid column")?;
    Ok((accession, taxid))
}

fn scan_rows(path: &Path, mut visit: impl FnMut(&str, TaxonId) -> Result<bool>) -> Result<()> {
    let mut reader = BufReader::new(GzDecoder::new(File::open(path)?));
    let (accession_column, taxid_column) = a2t_columns(&mut reader)?;
    for line in reader.lines() {
        let line = line?;
        let mut accession = None;
        let mut taxid = None;
        for (index, field) in line.split('\t').enumerate() {
            if index == accession_column {
                accession = Some(field);
            }
            if index == taxid_column {
                taxid = Some(field);
            }
            if accession.is_some() && taxid.is_some() {
                break;
            }
        }
        let (Some(accession), Some(taxid)) = (accession, taxid) else {
            continue;
        };
        if !visit(accession, taxid.parse()?)? {
            break;
        }
    }
    Ok(())
}
fn lookup_a2t(
    save_folder: &Path,
    accessions: &[impl AsRef<str> + Sync],
    low_memory: bool,
    wgs: bool,
) -> Result<HashMap<String, TaxonId>> {
    let requested = parse_accessions(accessions, true)
        .into_iter()
        .filter(|value| value != "NA")
        .collect::<HashSet<_>>();
    if low_memory {
        let partials = ensure_a2t_files(save_folder, wgs, false)?
            .par_iter()
            .map(|path| {
                let mut found = HashMap::new();
                scan_rows(path, |accession, taxid| {
                    if requested.contains(accession) {
                        found.insert(accession.to_owned(), taxid);
                    }
                    Ok(found.len() != requested.len())
                })?;
                Ok(found)
            })
            .collect::<Result<Vec<HashMap<String, TaxonId>>>>()?;
        return Ok(partials.into_iter().flatten().collect());
    }

    ensure_a2t_db(save_folder, false, wgs)?;
    let mut connection = Connection::open(save_folder.join(DB_FILE))?;
    connection.execute_batch(
        "DROP TABLE IF EXISTS temp.tmp_accs;
         CREATE TEMP TABLE tmp_accs (accession TEXT PRIMARY KEY);",
    )?;
    let transaction = connection.transaction()?;
    {
        let mut insert = transaction.prepare("INSERT OR IGNORE INTO tmp_accs VALUES (?1)")?;
        for accession in requested {
            insert.execute([accession])?;
        }
    }
    transaction.commit()?;
    let mut statement = connection.prepare(
        "SELECT t.accession, a.taxid
         FROM tmp_accs t JOIN a2t a ON t.accession = a.accession",
    )?;
    let rows = statement.query_map([], |row| {
        Ok((row.get::<_, String>(0)?, row.get::<_, TaxonId>(1)?))
    })?;
    rows.collect::<rusqlite::Result<HashMap<_, _>>>()
        .map_err(Into::into)
}
fn lookup_t2a(
    save_folder: &Path,
    taxa: &[TaxonId],
    low_memory: bool,
    wgs: bool,
) -> Result<HashSet<String>> {
    let requested = taxa.iter().copied().collect::<HashSet<_>>();
    if low_memory {
        let partials = ensure_a2t_files(save_folder, wgs, false)?
            .par_iter()
            .map(|path| {
                let mut found = HashSet::new();
                scan_rows(path, |accession, taxid| {
                    if requested.contains(&taxid) {
                        found.insert(accession.to_owned());
                    }
                    Ok(true)
                })?;
                Ok(found)
            })
            .collect::<Result<Vec<HashSet<String>>>>()?;
        return Ok(partials.into_iter().flatten().collect());
    }

    ensure_a2t_db(save_folder, false, wgs)?;
    let mut connection = Connection::open(save_folder.join(DB_FILE))?;
    connection.execute_batch(
        "DROP TABLE IF EXISTS temp.tmp_taxa;
         CREATE TEMP TABLE tmp_taxa (taxid INTEGER PRIMARY KEY);",
    )?;
    let transaction = connection.transaction()?;
    {
        let mut insert = transaction.prepare("INSERT OR IGNORE INTO tmp_taxa VALUES (?1)")?;
        for taxid in requested {
            insert.execute([taxid])?;
        }
    }
    transaction.commit()?;
    let mut statement = connection
        .prepare("SELECT accession FROM a2t JOIN tmp_taxa ON a2t.taxid = tmp_taxa.taxid")?;
    let rows = statement.query_map([], |row| row.get::<_, String>(0))?;
    rows.collect::<rusqlite::Result<HashSet<_>>>()
        .map_err(Into::into)
}

fn ensure_a2t_db(save_folder: &Path, rebuild: bool, wgs: bool) -> Result<()> {
    let sources = ensure_a2t_files(save_folder, wgs, rebuild)?;
    let db_path = save_folder.join(DB_FILE);
    if rebuild && db_path.exists() {
        fs::remove_file(&db_path)?;
    }
    let mut connection = Connection::open(&db_path)?;
    connection.execute_batch(
        "CREATE TABLE IF NOT EXISTS a2t (accession TEXT, taxid INTEGER);
         CREATE TABLE IF NOT EXISTS a2t_sources (source TEXT PRIMARY KEY, status TEXT);",
    )?;
    let incomplete: i64 = connection.query_row(
        "SELECT COUNT(*) FROM a2t_sources WHERE status != 'complete'",
        [],
        |row| row.get(0),
    )?;
    if incomplete > 0 {
        drop(connection);
        fs::remove_file(&db_path)?;
        return ensure_a2t_db(save_folder, false, wgs);
    }
    let mut loaded = HashSet::new();
    {
        let mut statement =
            connection.prepare("SELECT source FROM a2t_sources WHERE status = 'complete'")?;
        let rows = statement.query_map([], |row| row.get::<_, String>(0))?;
        for row in rows {
            loaded.insert(row?);
        }
    }
    let existing_rows: i64 =
        connection.query_row("SELECT COUNT(*) FROM a2t", [], |row| row.get(0))?;
    if loaded.is_empty() && existing_rows > 0 {
        connection.execute(
            "INSERT OR REPLACE INTO a2t_sources(source, status) VALUES (?1, 'complete')",
            [GB_FILE],
        )?;
        loaded.insert(GB_FILE.to_owned());
        if fs::metadata(&db_path)?.len() >= 10_000_000_000 {
            connection.execute(
                "INSERT OR REPLACE INTO a2t_sources(source, status) VALUES (?1, 'complete')",
                [WGS_FILE],
            )?;
            loaded.insert(WGS_FILE.to_owned());
        }
    }
    for source in sources {
        let filename = source.file_name().unwrap().to_string_lossy().into_owned();
        if loaded.contains(&filename) {
            continue;
        }
        connection.execute(
            "INSERT OR REPLACE INTO a2t_sources(source, status) VALUES (?1, 'loading')",
            [&filename],
        )?;
        let transaction = connection.transaction()?;
        {
            let mut insert = transaction.prepare("INSERT INTO a2t VALUES (?1, ?2)")?;
            scan_rows(&source, |accession, taxid| {
                insert.execute(params![accession, taxid])?;
                Ok(true)
            })?;
        }
        transaction.commit()?;
        connection.execute(
            "UPDATE a2t_sources SET status = 'complete' WHERE source = ?1",
            [&filename],
        )?;
    }
    connection.execute_batch(
        "CREATE INDEX IF NOT EXISTS idx_accession ON a2t(accession);
         CREATE INDEX IF NOT EXISTS idx_taxid ON a2t(taxid);",
    )?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_small_dump_and_corrects_ranks() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("nodes.dmp");
        fs::write(
            &path,
            "1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tsuperkingdom\t|\n3\t|\t2\t|\tno rank\t|\n",
        )
        .unwrap();
        let nodes = build_nodes(&path).unwrap();
        assert_eq!(nodes[0].rank_code, "R");
        assert_eq!(nodes[1].rank_code, "D");
        assert_eq!(nodes[2].rank_code, "D2");
    }

    #[test]
    fn full_fixture_matches_python_reference() {
        let dir = tempfile::tempdir().unwrap();
        fs::write(
            dir.path().join("names.dmp"),
            concat!(
                "1\t|\troot\t|\t\t|\tscientific name\t|\n",
                "2\t|\tBacteria\t|\t\t|\tscientific name\t|\n",
                "10\t|\tFam\t|\t\t|\tscientific name\t|\n",
                "11\t|\tGen\t|\t\t|\tscientific name\t|\n",
                "12\t|\tSp A\t|\t\t|\tscientific name\t|\n",
                "13\t|\tSp B\t|\t\t|\tscientific name\t|\n",
                "14\t|\tStrain\t|\t\t|\tscientific name\t|\n",
            ),
        )
        .unwrap();
        fs::write(
            dir.path().join("nodes.dmp"),
            concat!(
                "1\t|\t1\t|\tno rank\t|\n",
                "2\t|\t1\t|\tsuperkingdom\t|\n",
                "10\t|\t2\t|\tfamily\t|\n",
                "11\t|\t10\t|\tgenus\t|\n",
                "12\t|\t11\t|\tspecies\t|\n",
                "13\t|\t11\t|\tspecies\t|\n",
                "14\t|\t12\t|\tno rank\t|\n",
            ),
        )
        .unwrap();
        fs::write(dir.path().join("targets.json"), r#"{"pathogens":{"a":12}}"#).unwrap();
        let tu = TaxonomicUtils::new(TaxutilsOptions {
            save_folder: dir.path().to_owned(),
            ..Default::default()
        })
        .unwrap();

        assert_eq!(tu.node(14).unwrap().rank_code, "S2");
        assert_eq!(tu.target_taxa, vec![10, 11, 12, 14]);
        assert_eq!(tu.get_branch(14), vec![1, 2, 10, 11, 12, 14]);
        assert_eq!(tu.get_subtree(11), vec![11, 12, 14, 13]);
        assert_eq!(tu.get_lca(14, 13), 11);
        assert_eq!(tu.get_distance(14, 13), 3);
        assert_eq!(tu.sort_taxa([14, 10, 13, 12]), vec![10, 12, 14, 13]);
        let profile = tu.topology(11, None).unwrap();
        assert_eq!(
            (profile.n_taxa, profile.n_leaves, profile.max_depth),
            (4, 2, 2)
        );
        assert_eq!(profile.mean_depth, 1.0);
        assert_eq!((profile.topology_scale, profile.max_children), (1, 2));
        assert_eq!(profile.branching_taxa_fraction, 0.5);
        assert_eq!(profile.top_child_fraction, 2.0 / 3.0);
    }

    fn write_a2t_fixture(path: &Path, rows: &[(&str, i64)]) {
        let file = File::create(path).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
        writeln!(encoder, "accession\taccession.version\ttaxid\tgi").unwrap();
        for (accession, taxid) in rows {
            writeln!(
                encoder,
                "{}\t{}\t{}\t0",
                accession.split('.').next().unwrap(),
                accession,
                taxid
            )
            .unwrap();
        }
        encoder.finish().unwrap();
    }

    #[test]
    fn parallel_gzip_and_batched_sqlite_lookups_are_equivalent() {
        let dir = tempfile::tempdir().unwrap();
        write_a2t_fixture(
            &dir.path().join(GB_FILE),
            &[("NC_000001.1", 10), ("NC_000002.1", 20)],
        );
        write_a2t_fixture(
            &dir.path().join(WGS_FILE),
            &[("ABCD01000001.1", 30), ("ABCD01000002.1", 20)],
        );
        let accessions = vec![
            "NC_000001.1".to_owned(),
            "ABCD01000001.1".to_owned(),
            "missing".to_owned(),
        ];
        let scanned = lookup_a2t(dir.path(), &accessions, true, true).unwrap();
        let indexed = lookup_a2t(dir.path(), &accessions, false, true).unwrap();
        assert_eq!(scanned, indexed);
        assert_eq!(scanned["NC_000001.1"], 10);
        assert_eq!(scanned["ABCD01000001.1"], 30);

        let scanned_reverse = lookup_t2a(dir.path(), &[20, 30], true, true).unwrap();
        let indexed_reverse = lookup_t2a(dir.path(), &[20, 30], false, true).unwrap();
        assert_eq!(scanned_reverse, indexed_reverse);
        assert_eq!(scanned_reverse.len(), 3);
    }
}
