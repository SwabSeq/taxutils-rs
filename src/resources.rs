use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use flate2::read::{GzDecoder, MultiGzDecoder};
use rapidgzip_core::Decoder;
use rayon::prelude::*;
use reqwest::blocking::Client;
use rusqlite::{Connection, params};
use serde_json::Value;

use crate::accession::parse_accessions;
use crate::taxonomy::{
    AccessionLookupOptions, TaxonId, TaxonNode, TaxonomicUtils, assign_rank_codes, canonical_name,
    rank_index, taxonomic_order,
};

const TAXDUMP_URL: &str = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz";
const TARGETS_URL: &str = "https://web.cs.ucla.edu/~wob/projects/taxutils/targets.json";
const A2T_BASE_URL: &str = "https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid";
const GB_FILE: &str = "nucl_gb.accession2taxid.gz";
const WGS_FILE: &str = "nucl_wgs.accession2taxid.gz";
const DB_FILE: &str = "nucl.accession2taxid.db";
const A2T_JOIN_SQL: &str = "SELECT t.accession, a.taxid
     FROM tmp_accs t
     CROSS JOIN a2t a INDEXED BY idx_accession ON t.accession = a.accession";
// Drive reverse lookups from the requested taxa. An ordinary JOIN can make
// SQLite scan the entire accession table and probe tmp_taxa for every row.
const T2A_JOIN_SQL: &str = "SELECT a.accession
     FROM tmp_taxa t
     CROSS JOIN a2t a INDEXED BY idx_taxid ON a.taxid = t.taxid";

/// Controls creation and retention of the indexed accession database.
#[derive(Clone, Copy, Debug)]
pub struct AccessionDatabaseOptions {
    pub rebuild: bool,
    pub wgs: bool,
    /// Retain the compressed NCBI input files after the SQLite database is ready.
    pub keep_downloads: bool,
}

impl Default for AccessionDatabaseOptions {
    fn default() -> Self {
        Self {
            rebuild: false,
            wgs: false,
            keep_downloads: true,
        }
    }
}

#[derive(Clone, Debug)]
pub struct TaxutilsOptions {
    pub accessions: Option<Vec<String>>,
    pub low_memory: bool,
    pub targets_json: Option<PathBuf>,
    pub rebuild: bool,
    pub wgs: bool,
    pub keep_accession_downloads: bool,
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
            keep_accession_downloads: true,
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
    /// Keep the compressed NCBI accession sources after the SQLite index is ready.
    pub fn keep_accession_downloads(mut self, value: bool) -> Self {
        self.options.keep_accession_downloads = value;
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
        ensure_a2t_db(
            &options.save_folder,
            options.rebuild,
            options.wgs,
            options.keep_accession_downloads,
        )?;
    }
    let mut a2t = HashMap::new();
    if let Some(accessions) = &options.accessions {
        a2t = lookup_a2t(
            &options.save_folder,
            accessions,
            options.low_memory,
            options.wgs,
            options.keep_accession_downloads,
        )?;
    }
    names.insert(2697049, "SARS-CoV-2".to_owned());
    names.insert(694009, "SARS-related-CoV".to_owned());
    Ok(TaxonomicUtils::from_parts(
        names,
        nodes,
        target_taxa,
        a2t,
        AccessionLookupOptions {
            low_memory: options.low_memory,
            wgs: options.wgs,
            keep_downloads: options.keep_accession_downloads,
            save_folder: options.save_folder,
        },
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
            self.keep_accession_downloads,
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
            self.keep_accession_downloads,
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
    ensure_a2t_paths(&paths, rebuild)?;
    Ok(paths)
}

fn ensure_a2t_paths(paths: &[PathBuf], rebuild: bool) -> Result<()> {
    paths.par_iter().try_for_each(|path| {
        if rebuild || !path.exists() {
            let filename = path.file_name().unwrap().to_string_lossy();
            download_file(&format!("{A2T_BASE_URL}/{filename}"), path)?;
        }
        Ok::<_, anyhow::Error>(())
    })?;
    Ok(())
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
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(path)?));
    scan_reader(&mut reader, &mut visit)
}

// Decompression overlaps with row scanning on a separate thread. Wider
// speculative decoder pools were slower on NCBI mapping data and could use
// hundreds of MiB during startup. One decoder per source keeps the pipeline
// small; GB and WGS still run concurrently in the surrounding Rayon pool.
fn lookup_decoder() -> Result<Decoder> {
    Ok(Decoder::builder()
        .decoder_threads(1)
        .decoded_chunk_size(1 << 20)
        .in_flight_chunks(2)
        .build()?)
}

fn scan_lookup_rows(
    path: &Path,
    decoder: &Decoder,
    mut visit: impl FnMut(&str, TaxonId) -> Result<bool>,
) -> Result<()> {
    let mut reader = BufReader::with_capacity(1 << 20, decoder.open(path)?);
    scan_reader(&mut reader, &mut visit)
        .with_context(|| format!("failed to scan {}", path.display()))
}

fn scan_reader(
    reader: &mut impl BufRead,
    visit: &mut impl FnMut(&str, TaxonId) -> Result<bool>,
) -> Result<()> {
    let (accession_column, taxid_column) = a2t_columns(reader)?;
    let mut visit_line = |line: &str| -> Result<bool> {
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
            return Ok(true);
        };
        visit(accession, taxid.parse()?)
    };
    // Borrow complete rows directly from the read buffer. Only a row crossing
    // a buffer boundary needs copying, into a reusable scratch allocation.
    let mut partial = Vec::new();
    loop {
        let bytes = reader.fill_buf()?;
        if bytes.is_empty() {
            if !partial.is_empty() {
                visit_line(std::str::from_utf8(&partial)?.trim_end_matches('\r'))?;
            }
            return Ok(());
        }
        let mut start = 0;
        if !partial.is_empty() {
            if let Some(end) = memchr::memchr(b'\n', bytes) {
                partial.extend_from_slice(&bytes[..end]);
                if !visit_line(std::str::from_utf8(&partial)?.trim_end_matches('\r'))? {
                    return Ok(());
                }
                partial.clear();
                start = end + 1;
            } else {
                partial.extend_from_slice(bytes);
                let consumed = bytes.len();
                reader.consume(consumed);
                continue;
            }
        }
        if let Some(end) = memchr::memrchr(b'\n', &bytes[start..]) {
            let end = start + end + 1;
            for line in std::str::from_utf8(&bytes[start..end])?.lines() {
                if !visit_line(line)? {
                    return Ok(());
                }
            }
            start = end;
        }
        partial.extend_from_slice(&bytes[start..]);
        let consumed = bytes.len();
        reader.consume(consumed);
    }
}
fn lookup_a2t(
    save_folder: &Path,
    accessions: &[impl AsRef<str> + Sync],
    low_memory: bool,
    wgs: bool,
    keep_accession_downloads: bool,
) -> Result<HashMap<String, TaxonId>> {
    let requested = parse_accessions(accessions, true)
        .into_iter()
        .filter(|value| value != "NA")
        .collect::<HashSet<_>>();
    lookup_parsed_a2t(
        save_folder,
        requested,
        low_memory,
        wgs,
        keep_accession_downloads,
    )
}

fn lookup_parsed_a2t(
    save_folder: &Path,
    requested: HashSet<String>,
    low_memory: bool,
    wgs: bool,
    keep_accession_downloads: bool,
) -> Result<HashMap<String, TaxonId>> {
    if requested.is_empty() {
        return Ok(HashMap::new());
    }
    if low_memory {
        let paths = ensure_a2t_files(save_folder, wgs, false)?;
        let decoder = lookup_decoder()?;
        let partials = paths
            .par_iter()
            .map(|path| {
                let mut found = HashMap::new();
                scan_lookup_rows(path, &decoder, |accession, taxid| {
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

    let mut index = AccessionTaxidIndex::open(save_folder, wgs, keep_accession_downloads)?;
    index.lookup(requested)
}

/// Reusable bounded-memory SQLite lookup handle for streaming FASTA workflows.
pub(crate) struct AccessionTaxidIndex {
    connection: Connection,
}

impl AccessionTaxidIndex {
    pub(crate) fn open(
        save_folder: impl AsRef<Path>,
        wgs: bool,
        keep_accession_downloads: bool,
    ) -> Result<Self> {
        let save_folder = save_folder.as_ref();
        fs::create_dir_all(save_folder)?;
        ensure_a2t_db(save_folder, false, wgs, keep_accession_downloads)?;
        Ok(Self {
            connection: Connection::open(save_folder.join(DB_FILE))?,
        })
    }

    pub(crate) fn lookup(
        &mut self,
        requested: HashSet<String>,
    ) -> Result<HashMap<String, TaxonId>> {
        lookup_a2t_connection(&mut self.connection, requested)
    }
}

fn lookup_a2t_connection(
    connection: &mut Connection,
    requested: HashSet<String>,
) -> Result<HashMap<String, TaxonId>> {
    connection.execute_batch(
        "PRAGMA temp_store = MEMORY;
         DROP TABLE IF EXISTS temp.tmp_accs;
         CREATE TEMP TABLE tmp_accs (accession TEXT PRIMARY KEY) WITHOUT ROWID;",
    )?;
    let transaction = connection.transaction()?;
    {
        let mut insert = transaction.prepare("INSERT INTO tmp_accs VALUES (?1)")?;
        for accession in requested {
            insert.execute([accession])?;
        }
    }
    transaction.commit()?;
    let mut statement = connection.prepare(A2T_JOIN_SQL)?;
    let rows = statement.query_map([], |row| {
        Ok((row.get::<_, String>(0)?, row.get::<_, TaxonId>(1)?))
    })?;
    rows.collect::<rusqlite::Result<HashMap<_, _>>>()
        .map_err(Into::into)
}

/// Look up parsed accession-to-taxid mappings without loading the taxonomy tree.
///
/// The accessions must already be normalized, versioned NCBI accessions. This is
/// useful for bulk FASTA workflows that only need the accession index and have
/// already parsed and deduplicated their headers.
pub fn lookup_accession_taxids(
    save_folder: impl AsRef<Path>,
    accessions: HashSet<String>,
    low_memory: bool,
    wgs: bool,
) -> Result<HashMap<String, TaxonId>> {
    let save_folder = save_folder.as_ref();
    fs::create_dir_all(save_folder)?;
    lookup_parsed_a2t(save_folder, accessions, low_memory, wgs, true)
}

/// Look up accessions assigned directly to the requested taxids without loading
/// the taxonomy tree.
pub fn lookup_taxid_accessions(
    save_folder: impl AsRef<Path>,
    taxa: &[TaxonId],
    low_memory: bool,
    wgs: bool,
) -> Result<HashSet<String>> {
    let save_folder = save_folder.as_ref();
    fs::create_dir_all(save_folder)?;
    lookup_t2a(save_folder, taxa, low_memory, wgs, true)
}

fn lookup_t2a(
    save_folder: &Path,
    taxa: &[TaxonId],
    low_memory: bool,
    wgs: bool,
    keep_accession_downloads: bool,
) -> Result<HashSet<String>> {
    let requested = taxa.iter().copied().collect::<HashSet<_>>();
    if requested.is_empty() {
        return Ok(HashSet::new());
    }
    if low_memory {
        let paths = ensure_a2t_files(save_folder, wgs, false)?;
        let decoder = lookup_decoder()?;
        let partials = paths
            .par_iter()
            .map(|path| {
                let mut found = HashSet::new();
                scan_lookup_rows(path, &decoder, |accession, taxid| {
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

    ensure_a2t_db(save_folder, false, wgs, keep_accession_downloads)?;
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
    let mut statement = connection.prepare(T2A_JOIN_SQL)?;
    let rows = statement.query_map([], |row| row.get::<_, String>(0))?;
    rows.collect::<rusqlite::Result<HashSet<_>>>()
        .map_err(Into::into)
}

#[derive(Debug)]
struct AccessionDatabaseState {
    loaded_sources: HashSet<String>,
    rebuild: bool,
}

fn requested_a2t_sources(wgs: bool) -> Vec<&'static str> {
    let mut sources = vec![GB_FILE];
    if wgs {
        sources.push(WGS_FILE);
    }
    sources
}

fn inspect_a2t_database(db_path: &Path) -> Result<AccessionDatabaseState> {
    let connection = Connection::open(db_path)?;
    connection.execute_batch(
        "CREATE TABLE IF NOT EXISTS a2t (accession TEXT, taxid INTEGER);
         CREATE TABLE IF NOT EXISTS a2t_sources (source TEXT PRIMARY KEY, status TEXT);",
    )?;
    let incomplete: bool = connection.query_row(
        "SELECT EXISTS(SELECT 1 FROM a2t_sources WHERE status != 'complete')",
        [],
        |row| row.get(0),
    )?;
    if incomplete {
        return Ok(AccessionDatabaseState {
            loaded_sources: HashSet::new(),
            rebuild: true,
        });
    }

    let has_rows: bool =
        connection.query_row("SELECT EXISTS(SELECT 1 FROM a2t LIMIT 1)", [], |row| {
            row.get(0)
        })?;
    if !has_rows {
        return Ok(AccessionDatabaseState {
            loaded_sources: HashSet::new(),
            rebuild: true,
        });
    }

    let mut loaded_sources = HashSet::new();
    {
        let mut statement =
            connection.prepare("SELECT source FROM a2t_sources WHERE status = 'complete'")?;
        let rows = statement.query_map([], |row| row.get::<_, String>(0))?;
        for row in rows {
            loaded_sources.insert(row?);
        }
    }
    if loaded_sources.is_empty() {
        connection.execute(
            "INSERT OR REPLACE INTO a2t_sources(source, status) VALUES (?1, 'complete')",
            [GB_FILE],
        )?;
        loaded_sources.insert(GB_FILE.to_owned());
        if fs::metadata(db_path)?.len() >= 10_000_000_000 {
            connection.execute(
                "INSERT OR REPLACE INTO a2t_sources(source, status) VALUES (?1, 'complete')",
                [WGS_FILE],
            )?;
            loaded_sources.insert(WGS_FILE.to_owned());
        }
    }
    connection.execute_batch(
        "CREATE INDEX IF NOT EXISTS idx_accession ON a2t(accession);
         CREATE INDEX IF NOT EXISTS idx_taxid ON a2t(taxid);",
    )?;
    Ok(AccessionDatabaseState {
        loaded_sources,
        rebuild: false,
    })
}

fn configure_bulk_load(connection: &Connection) -> Result<()> {
    let worker_threads = std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
        .min(8);
    connection.pragma_update(None, "journal_mode", "OFF")?;
    connection.pragma_update(None, "synchronous", "OFF")?;
    connection.pragma_update(None, "locking_mode", "EXCLUSIVE")?;
    connection.pragma_update(None, "cache_size", -131_072_i64)?;
    connection.pragma_update(None, "threads", worker_threads as i64)?;
    Ok(())
}

fn insert_a2t_source(connection: &mut Connection, source: &Path) -> Result<()> {
    let filename = source.file_name().unwrap().to_string_lossy().into_owned();
    connection.execute(
        "INSERT OR REPLACE INTO a2t_sources(source, status) VALUES (?1, 'loading')",
        [&filename],
    )?;
    let transaction = connection.transaction()?;
    {
        let mut insert = transaction.prepare("INSERT INTO a2t VALUES (?1, ?2)")?;
        scan_rows(source, |accession, taxid| {
            insert.execute(params![accession, taxid])?;
            Ok(true)
        })?;
    }
    transaction.commit()?;
    connection.execute(
        "UPDATE a2t_sources SET status = 'complete' WHERE source = ?1",
        [&filename],
    )?;
    Ok(())
}

fn finish_bulk_load(connection: &Connection) -> Result<()> {
    connection.execute_batch(
        "CREATE INDEX idx_accession ON a2t(accession);
         CREATE INDEX idx_taxid ON a2t(taxid);",
    )?;
    connection.pragma_update(None, "locking_mode", "NORMAL")?;
    connection.pragma_update(None, "journal_mode", "DELETE")?;
    connection.pragma_update(None, "synchronous", "FULL")?;
    Ok(())
}

fn install_temporary_database(temporary: tempfile::TempPath, db_path: &Path) -> Result<()> {
    File::open(&temporary)?.sync_all()?;
    temporary
        .persist(db_path)
        .map_err(|error| error.error)
        .with_context(|| format!("failed to install {}", db_path.display()))?;
    #[cfg(unix)]
    File::open(db_path.parent().unwrap_or_else(|| Path::new(".")))?.sync_all()?;
    Ok(())
}

fn build_a2t_database_atomic(sources: &[PathBuf], db_path: &Path) -> Result<()> {
    let parent = db_path.parent().unwrap_or_else(|| Path::new("."));
    let temporary = tempfile::NamedTempFile::new_in(parent)?.into_temp_path();
    let mut connection = Connection::open(&temporary)?;
    configure_bulk_load(&connection)?;
    connection.execute_batch(
        "CREATE TABLE a2t (accession TEXT, taxid INTEGER);
         CREATE TABLE a2t_sources (source TEXT PRIMARY KEY, status TEXT);",
    )?;
    for source in sources {
        insert_a2t_source(&mut connection, source)?;
    }
    finish_bulk_load(&connection)?;
    connection.close().map_err(|(_, error)| error)?;
    install_temporary_database(temporary, db_path)
}

fn upgrade_a2t_database_atomic(db_path: &Path, missing_sources: &[PathBuf]) -> Result<()> {
    let parent = db_path.parent().unwrap_or_else(|| Path::new("."));
    let temporary = tempfile::NamedTempFile::new_in(parent)?.into_temp_path();
    fs::copy(db_path, &temporary)?;
    let mut connection = Connection::open(&temporary)?;
    configure_bulk_load(&connection)?;
    connection.execute_batch(
        "DROP INDEX IF EXISTS idx_accession;
         DROP INDEX IF EXISTS idx_taxid;",
    )?;
    for source in missing_sources {
        insert_a2t_source(&mut connection, source)?;
    }
    finish_bulk_load(&connection)?;
    connection.close().map_err(|(_, error)| error)?;
    install_temporary_database(temporary, db_path)
}

fn discard_a2t_downloads(save_folder: &Path, wgs: bool) -> Result<()> {
    for path in a2t_paths(save_folder, wgs) {
        if path.exists() {
            fs::remove_file(&path)
                .with_context(|| format!("failed to remove {}", path.display()))?;
        }
    }
    Ok(())
}

fn ensure_a2t_db(save_folder: &Path, rebuild: bool, wgs: bool, keep_downloads: bool) -> Result<()> {
    fs::create_dir_all(save_folder)?;
    let db_path = save_folder.join(DB_FILE);
    let requested_sources = requested_a2t_sources(wgs);

    if !rebuild && db_path.exists() {
        let state = inspect_a2t_database(&db_path)?;
        if !state.rebuild {
            let missing_paths = requested_sources
                .iter()
                .filter(|source| !state.loaded_sources.contains(**source))
                .map(|source| save_folder.join(source))
                .collect::<Vec<_>>();
            if missing_paths.is_empty() {
                if !keep_downloads {
                    discard_a2t_downloads(save_folder, wgs)?;
                }
                return Ok(());
            }
            ensure_a2t_paths(&missing_paths, false)?;
            upgrade_a2t_database_atomic(&db_path, &missing_paths)?;
            if !keep_downloads {
                discard_a2t_downloads(save_folder, wgs)?;
            }
            return Ok(());
        }
    }

    let sources = ensure_a2t_files(save_folder, wgs, rebuild)?;
    build_a2t_database_atomic(&sources, &db_path)?;
    if !keep_downloads {
        discard_a2t_downloads(save_folder, wgs)?;
    }
    Ok(())
}

/// Build or validate the shared SQLite accession index.
///
/// New databases and upgrades are assembled beside the destination and
/// atomically installed only after all rows and indexes are complete.
pub fn ensure_accession_database(
    save_folder: impl AsRef<Path>,
    options: AccessionDatabaseOptions,
) -> Result<PathBuf> {
    let save_folder = save_folder.as_ref();
    ensure_a2t_db(
        save_folder,
        options.rebuild,
        options.wgs,
        options.keep_downloads,
    )?;
    Ok(save_folder.join(DB_FILE))
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
    fn buffered_scan_preserves_rows_at_every_boundary() {
        let input = b"taxid\tgi\taccession.version\taccession\r\n13\t0\tNC_000001.1\tNC_000001\r\nshort\n15\t0\tNC_000002.1\tNC_000002";
        for capacity in 1..=input.len() {
            let mut reader = BufReader::with_capacity(capacity, &input[..]);
            let mut rows = Vec::new();
            scan_reader(&mut reader, &mut |accession, taxid| {
                rows.push((accession.to_owned(), taxid));
                Ok(true)
            })
            .unwrap();
            assert_eq!(
                rows,
                vec![
                    ("NC_000001.1".to_owned(), 13),
                    ("NC_000002.1".to_owned(), 15)
                ],
                "capacity {capacity}"
            );
        }
    }

    #[test]
    fn empty_lookups_do_not_download_or_build() {
        let dir = tempfile::tempdir().unwrap();
        for low_memory in [true, false] {
            assert!(
                lookup_accession_taxids(dir.path(), HashSet::new(), low_memory, true)
                    .unwrap()
                    .is_empty()
            );
            assert!(
                lookup_taxid_accessions(dir.path(), &[], low_memory, true)
                    .unwrap()
                    .is_empty()
            );
        }
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
    }

    #[test]
    fn lookup_reads_concatenated_members_and_honors_wgs() {
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join(GB_FILE);
        write_a2t_fixture(&source, &[("NC_000001.1", 13)]);
        let file = fs::OpenOptions::new().append(true).open(&source).unwrap();
        let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::fast());
        write!(encoder, "NC_000002\tNC_000002.1\t15\t0").unwrap();
        encoder.finish().unwrap();
        write_a2t_fixture(
            &dir.path().join(WGS_FILE),
            &[("ABCD01000001.1", 15), ("NC_000001.1", 99)],
        );
        let queries = HashSet::from([
            "NC_000001.1".to_owned(),
            "NC_000002.1".to_owned(),
            "ABCD01000001.1".to_owned(),
        ]);
        for threads in [1, 4] {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    let gb =
                        lookup_accession_taxids(dir.path(), queries.clone(), true, false).unwrap();
                    assert_eq!(
                        gb,
                        HashMap::from([
                            ("NC_000001.1".to_owned(), 13),
                            ("NC_000002.1".to_owned(), 15)
                        ])
                    );
                    let both =
                        lookup_accession_taxids(dir.path(), queries.clone(), true, true).unwrap();
                    assert_eq!(both["NC_000001.1"], 99); // WGS keeps its existing precedence.
                    assert_eq!(both["ABCD01000001.1"], 15);
                    assert_eq!(
                        lookup_taxid_accessions(dir.path(), &[15, 15], true, false).unwrap(),
                        HashSet::from(["NC_000002.1".to_owned()])
                    );
                    assert_eq!(
                        lookup_taxid_accessions(dir.path(), &[15], true, true).unwrap(),
                        HashSet::from(["NC_000002.1".to_owned(), "ABCD01000001.1".to_owned()])
                    );
                });
        }
    }

    #[test]
    fn full_scans_reject_corrupt_or_truncated_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join(GB_FILE);
        write_a2t_fixture(&source, &[("NC_000001.1", 13)]);
        let original = fs::read(&source).unwrap();
        let mut corrupt = original.clone();
        let checksum = corrupt.len() - 8;
        corrupt[checksum] ^= 0xff;
        for invalid in [corrupt, original[..original.len() - 4].to_vec()] {
            fs::write(&source, invalid).unwrap();
            assert!(
                lookup_accession_taxids(
                    dir.path(),
                    HashSet::from(["missing".to_owned()]),
                    true,
                    false
                )
                .is_err()
            );
            assert!(lookup_taxid_accessions(dir.path(), &[13], true, false).is_err());
        }
    }

    #[test]
    fn parallel_scan_handles_multiple_buffers_and_early_stop() {
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join(GB_FILE);
        let mut encoder = flate2::write::GzEncoder::new(
            File::create(&source).unwrap(),
            flate2::Compression::fast(),
        );
        writeln!(encoder, "accession\taccession.version\ttaxid\tgi").unwrap();
        for index in 0..100_000 {
            writeln!(
                encoder,
                "NC_{index:09}\tNC_{index:09}.1\t{}\t0",
                index % 1000
            )
            .unwrap();
        }
        // Early stop must not visit rows beyond the complete accession result.
        writeln!(encoder, "NC_bad\tNC_bad.1\tinvalid_taxid\t0").unwrap();
        encoder.finish().unwrap();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap();
        pool.install(|| {
            let found = lookup_accession_taxids(
                dir.path(),
                HashSet::from(["NC_000000000.1".to_owned(), "NC_000099999.1".to_owned()]),
                true,
                false,
            )
            .unwrap();
            assert_eq!(found["NC_000000000.1"], 0);
            assert_eq!(found["NC_000099999.1"], 999);
            assert!(lookup_taxid_accessions(dir.path(), &[999], true, false).is_err());
        });
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
        let scanned = lookup_a2t(dir.path(), &accessions, true, true, true).unwrap();
        let indexed = lookup_a2t(dir.path(), &accessions, false, true, true).unwrap();
        let direct = lookup_accession_taxids(
            dir.path(),
            accessions.iter().cloned().collect(),
            false,
            true,
        )
        .unwrap();
        assert_eq!(scanned, indexed);
        assert_eq!(direct, indexed);
        assert_eq!(scanned["NC_000001.1"], 10);
        assert_eq!(scanned["ABCD01000001.1"], 30);

        let connection = Connection::open(dir.path().join(DB_FILE)).unwrap();
        connection
            .execute_batch("CREATE TEMP TABLE tmp_accs (accession TEXT PRIMARY KEY) WITHOUT ROWID;")
            .unwrap();
        let mut plan = connection
            .prepare(&format!("EXPLAIN QUERY PLAN {A2T_JOIN_SQL}"))
            .unwrap();
        let details = plan
            .query_map([], |row| row.get::<_, String>(3))
            .unwrap()
            .collect::<rusqlite::Result<Vec<_>>>()
            .unwrap();
        assert!(
            details
                .iter()
                .any(|detail| detail.contains("idx_accession")),
            "query plan did not use idx_accession: {details:?}"
        );

        let scanned_reverse = lookup_t2a(dir.path(), &[20, 30], true, true, true).unwrap();
        let indexed_reverse = lookup_t2a(dir.path(), &[20, 30], false, true, true).unwrap();
        assert_eq!(scanned_reverse, indexed_reverse);
        assert_eq!(scanned_reverse.len(), 3);
    }

    #[test]
    fn reverse_lookup_searches_taxid_index_and_preserves_set_semantics() {
        let dir = tempfile::tempdir().unwrap();
        write_a2t_fixture(
            &dir.path().join(GB_FILE),
            &[
                ("NC_000001.1", 12059),
                ("NC_000002.1", 12059),
                ("NC_000003.1", 20),
            ],
        );
        ensure_accession_database(dir.path(), AccessionDatabaseOptions::default()).unwrap();
        let connection = Connection::open(dir.path().join(DB_FILE)).unwrap();
        connection
            .execute_batch(
                "CREATE TEMP TABLE tmp_taxa (taxid INTEGER PRIMARY KEY);
                            INSERT INTO tmp_taxa VALUES (12059);",
            )
            .unwrap();
        let mut plan = connection
            .prepare(&format!("EXPLAIN QUERY PLAN {T2A_JOIN_SQL}"))
            .unwrap();
        let details = plan
            .query_map([], |row| row.get::<_, String>(3))
            .unwrap()
            .collect::<rusqlite::Result<Vec<_>>>()
            .unwrap();
        assert!(
            details.iter().any(|detail| {
                detail.contains("SEARCH a USING INDEX idx_taxid") && detail.contains("taxid=?")
            }),
            "reverse lookup must search idx_taxid, not scan a2t: {details:?}"
        );
        assert!(
            !details
                .iter()
                .any(|detail| detail == "SCAN a" || detail.starts_with("SCAN a ")),
            "reverse lookup scanned the accession table: {details:?}"
        );
        for taxa in [
            vec![12059],
            vec![12059, 12059, 99999],
            vec![12059, 20],
            vec![99999],
            vec![],
        ] {
            let scanned = lookup_taxid_accessions(dir.path(), &taxa, true, false).unwrap();
            let indexed = lookup_taxid_accessions(dir.path(), &taxa, false, false).unwrap();
            assert_eq!(indexed, scanned, "reverse lookup mismatch for {taxa:?}");
        }
        assert_eq!(
            lookup_taxid_accessions(dir.path(), &[12059], false, false).unwrap(),
            HashSet::from(["NC_000001.1".to_owned(), "NC_000002.1".to_owned()])
        );
    }

    #[test]
    fn database_build_is_atomic_and_can_discard_downloads() {
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join(GB_FILE);
        write_a2t_fixture(&source, &[("NC_000001.1", 10), ("NC_000002.1", 20)]);

        let db_path = ensure_accession_database(
            dir.path(),
            AccessionDatabaseOptions {
                keep_downloads: false,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(db_path.exists());
        assert!(!source.exists());
        assert!(
            fs::read_dir(dir.path()).unwrap().all(|entry| !entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .contains("tmp"))
        );

        // A complete database can be reused without downloading the discarded input.
        ensure_accession_database(
            dir.path(),
            AccessionDatabaseOptions {
                keep_downloads: false,
                ..Default::default()
            },
        )
        .unwrap();
        let connection = Connection::open(db_path).unwrap();
        let rows: i64 = connection
            .query_row("SELECT COUNT(*) FROM a2t", [], |row| row.get(0))
            .unwrap();
        assert_eq!(rows, 2);
        let indexes = connection
            .prepare("SELECT name FROM sqlite_master WHERE type = 'index' AND name LIKE 'idx_%'")
            .unwrap()
            .query_map([], |row| row.get::<_, String>(0))
            .unwrap()
            .collect::<rusqlite::Result<HashSet<_>>>()
            .unwrap();
        assert_eq!(
            indexes,
            HashSet::from(["idx_accession".to_owned(), "idx_taxid".to_owned()])
        );
    }

    #[test]
    fn failed_atomic_build_preserves_existing_database() {
        let dir = tempfile::tempdir().unwrap();
        let db_path = dir.path().join(DB_FILE);
        fs::write(&db_path, b"existing database sentinel").unwrap();
        let invalid_source = dir.path().join(GB_FILE);
        fs::write(
            &invalid_source,
            "accession\taccession.version\ttaxid\tgi\nNC_1\tNC_1.1\tnot-a-taxid\t0\n",
        )
        .unwrap();

        assert!(build_a2t_database_atomic(&[invalid_source], &db_path).is_err());
        assert_eq!(fs::read(&db_path).unwrap(), b"existing database sentinel");
    }

    #[test]
    fn empty_database_with_complete_metadata_is_rebuilt() {
        let dir = tempfile::tempdir().unwrap();
        let db_path = dir.path().join(DB_FILE);
        let connection = Connection::open(&db_path).unwrap();
        connection
            .execute_batch(
                "CREATE TABLE a2t (accession TEXT, taxid INTEGER);
                 CREATE TABLE a2t_sources (source TEXT PRIMARY KEY, status TEXT);
                 INSERT INTO a2t_sources VALUES ('nucl_gb.accession2taxid.gz', 'complete');",
            )
            .unwrap();
        drop(connection);
        write_a2t_fixture(&dir.path().join(GB_FILE), &[("NC_000001.1", 10)]);

        ensure_accession_database(dir.path(), AccessionDatabaseOptions::default()).unwrap();
        let connection = Connection::open(db_path).unwrap();
        let rows: i64 = connection
            .query_row("SELECT COUNT(*) FROM a2t", [], |row| row.get(0))
            .unwrap();
        assert_eq!(rows, 1);
    }
}
