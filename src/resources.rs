use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::thread;
use std::time::Duration;

use anyhow::{Context, Result, bail};
use flate2::read::GzDecoder;
use rapidgzip_core::Decoder;
use rayon::prelude::*;
use reqwest::blocking::Client;
use rusqlite::{Connection, params};
use serde_json::Value;

use crate::accession::parse_accessions;
use crate::fasta::CancellationToken;
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
/// Bumped whenever the on-disk layout changes in a way older builds cannot read.
/// v2 keys `a2t` on the accession itself (`WITHOUT ROWID`), dropping the separate
/// `idx_accession` and every stored rowid.
const SCHEMA_VERSION: i64 = 2;
/// Rows moved across the reader/merger channel at a time.
const MERGE_BATCH_ROWS: usize = 64 * 1024;
/// Rows scanned between cancellation checks, to keep the atomic load off the hot path.
const CANCEL_CHECK_ROWS: u64 = 64 * 1024;
const A2T_SCHEMA_SQL: &str = "\
     CREATE TABLE a2t (accession TEXT PRIMARY KEY, taxid INTEGER) WITHOUT ROWID;
     CREATE TABLE a2t_meta (key TEXT PRIMARY KEY, value TEXT);
     CREATE TABLE a2t_sources (
         source TEXT PRIMARY KEY,
         status TEXT,
         etag TEXT,
         last_modified TEXT,
         size INTEGER,
         row_count INTEGER
     );";
const A2T_JOIN_SQL: &str = "SELECT t.accession, a.taxid
     FROM tmp_accs t
     CROSS JOIN a2t a ON t.accession = a.accession";
// Drive reverse lookups from the requested taxa. An ordinary JOIN can make
// SQLite scan the entire accession table and probe tmp_taxa for every row.
const T2A_JOIN_SQL: &str = "SELECT a.accession
     FROM tmp_taxa t
     CROSS JOIN a2t a INDEXED BY idx_taxid ON a.taxid = t.taxid";

/// Controls creation and retention of the indexed accession database.
#[derive(Clone, Copy, Debug)]
pub struct AccessionDatabaseOptions {
    /// Ask NCBI whether the sources have changed and, if so, apply the
    /// difference to the existing database row by row.
    ///
    /// Off by default: opening a database must not depend on the network.
    pub refresh: bool,
    pub wgs: bool,
    /// Worker threads for every parallel stage; `None` means all logical CPUs.
    pub threads: Option<usize>,
}

impl Default for AccessionDatabaseOptions {
    fn default() -> Self {
        Self {
            refresh: false,
            wgs: false,
            threads: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct TaxutilsOptions {
    pub accessions: Option<Vec<String>>,
    pub low_memory: bool,
    pub targets_json: Option<PathBuf>,
    pub refresh: bool,
    pub wgs: bool,
    pub save_folder: PathBuf,
    /// Worker threads for every parallel stage; `None` means all logical CPUs.
    pub threads: Option<usize>,
}

impl Default for TaxutilsOptions {
    fn default() -> Self {
        Self {
            accessions: None,
            low_memory: true,
            targets_json: None,
            refresh: false,
            wgs: false,
            save_folder: std::env::var_os("TAXUTILS_GLOBALS")
                .map(PathBuf::from)
                .unwrap_or_else(|| PathBuf::from("./taxutils/")),
            threads: None,
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
    /// Re-fetch the managed taxonomy files and bring the accession database
    /// up to date, applying only the rows that changed upstream.
    pub fn refresh(mut self, value: bool) -> Self {
        self.options.refresh = value;
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
    /// Worker threads for every parallel stage. `None` means all logical CPUs.
    pub fn threads(mut self, value: Option<usize>) -> Self {
        self.options.threads = value;
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
    if options.refresh || !names_path.exists() || !nodes_path.exists() {
        download_taxonomy(&options.save_folder, &names_path, &nodes_path)?;
    }
    let targets_path = if let Some(path) = &options.targets_json {
        path.clone()
    } else {
        let path = options.save_folder.join("targets.json");
        if options.refresh || !path.exists() {
            download_file(TARGETS_URL, &path)?;
        }
        path
    };
    let mut names = build_names(&names_path)?;
    let nodes = build_nodes(&nodes_path)?;
    let target_taxa = build_target_taxa(&nodes, &names, &targets_path)?;
    if options.refresh || !options.low_memory {
        ensure_a2t_db(
            &options.save_folder,
            options.refresh,
            options.wgs,
            crate::threads::resolve(options.threads)?,
            &CancellationToken::default(),
        )?;
    }
    let mut a2t = HashMap::new();
    if let Some(accessions) = &options.accessions {
        a2t = lookup_a2t(
            &options.save_folder,
            accessions,
            options.low_memory,
            options.wgs,
            options.threads,
            &CancellationToken::default(),
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
            save_folder: options.save_folder,
            threads: options.threads,
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
            self.threads,
            &CancellationToken::default(),
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
            self.threads,
            &CancellationToken::default(),
        )
    }
}

fn download_file(url: &str, path: &Path) -> Result<()> {
    download_file_with_meta(url, path).map(|_| ())
}

fn download_file_with_meta(url: &str, path: &Path) -> Result<SourceMeta> {
    let parent = path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent)?;
    let mut temporary = tempfile::NamedTempFile::new_in(parent)?;
    let mut response = Client::new()
        .get(url)
        .send()
        .with_context(|| format!("failed to download {url}"))?
        .error_for_status()?;
    let meta = response_meta(&response);
    std::io::copy(&mut response, &mut temporary)?;
    temporary.flush()?;
    temporary
        .persist(path)
        .map_err(|error| error.error)
        .with_context(|| format!("failed to install {}", path.display()))?;
    Ok(meta)
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
fn ensure_a2t_files(save_folder: &Path, wgs: bool) -> Result<Vec<PathBuf>> {
    let paths = a2t_paths(save_folder, wgs);
    ensure_a2t_paths(&paths)?;
    Ok(paths)
}

fn ensure_a2t_paths(paths: &[PathBuf]) -> Result<()> {
    paths.par_iter().try_for_each(|path| {
        if !path.exists() {
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

// Decompression overlaps with row scanning on a separate thread. Wider
// speculative decoder pools were slower on NCBI mapping data and could use
// hundreds of MiB during startup. One decoder per source keeps the pipeline
// small; GB and WGS still run concurrently in the surrounding Rayon pool.
fn lookup_decoder(threads: usize) -> Result<Decoder> {
    Ok(Decoder::builder()
        .decoder_threads(threads.max(1))
        .decoded_chunk_size(1 << 20)
        .in_flight_chunks(2)
        .build()?)
}

// The database build is one pass over both sources and decompression is the
// limiting stage, so it uses larger chunks than `lookup_decoder`. The worker
// count comes from the caller either way.
fn build_decoder(threads: usize) -> Result<Decoder> {
    Ok(Decoder::builder()
        .decoder_threads(threads.max(1))
        .decoded_chunk_size(4 << 20)
        .in_flight_chunks(4)
        .build()?)
}

fn scan_lookup_rows(
    path: &Path,
    decoder: &Decoder,
    cancellation: &CancellationToken,
    mut visit: impl FnMut(&str, TaxonId) -> Result<bool>,
) -> Result<()> {
    let mut reader = BufReader::with_capacity(1 << 20, decoder.open(path)?);
    scan_reader(&mut reader, cancellation, &mut visit)
        .with_context(|| format!("failed to scan {}", path.display()))
}

fn scan_reader(
    reader: &mut impl BufRead,
    cancellation: &CancellationToken,
    visit: &mut impl FnMut(&str, TaxonId) -> Result<bool>,
) -> Result<()> {
    let (accession_column, taxid_column) = a2t_columns(reader)?;
    let mut scanned: u64 = 0;
    let mut visit_line = |line: &str| -> Result<bool> {
        scanned += 1;
        if scanned.is_multiple_of(CANCEL_CHECK_ROWS) && cancellation.is_cancelled() {
            bail!("operation cancelled");
        }
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
    threads: Option<usize>,
    cancellation: &CancellationToken,
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
        threads,
        cancellation,
    )
}

fn lookup_parsed_a2t(
    save_folder: &Path,
    requested: HashSet<String>,
    low_memory: bool,
    wgs: bool,
    threads: Option<usize>,
    cancellation: &CancellationToken,
) -> Result<HashMap<String, TaxonId>> {
    if requested.is_empty() {
        return Ok(HashMap::new());
    }
    if low_memory {
        let paths = ensure_a2t_files(save_folder, wgs)?;
        let decoder = lookup_decoder(crate::threads::resolve(threads)?)?;
        let partials = paths
            .par_iter()
            .map(|path| {
                let mut found = HashMap::new();
                scan_lookup_rows(path, &decoder, cancellation, |accession, taxid| {
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

    let mut index =
        AccessionTaxidIndex::open(save_folder, wgs, threads, cancellation)?;
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
        threads: Option<usize>,
        cancellation: &CancellationToken,
    ) -> Result<Self> {
        let save_folder = save_folder.as_ref();
        fs::create_dir_all(save_folder)?;
        ensure_a2t_db(
            save_folder,
            false,
            wgs,
            crate::threads::resolve(threads)?,
            cancellation,
        )?;
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
    threads: Option<usize>,
) -> Result<HashMap<String, TaxonId>> {
    lookup_accession_taxids_with_cancel(
        save_folder,
        accessions,
        low_memory,
        wgs,
        threads,
        &CancellationToken::default(),
    )
}

/// [`lookup_accession_taxids`] that stops early when `cancellation` is triggered.
///
/// A build started on behalf of this lookup is abandoned without touching any
/// database already installed at the destination.
pub fn lookup_accession_taxids_with_cancel(
    save_folder: impl AsRef<Path>,
    accessions: HashSet<String>,
    low_memory: bool,
    wgs: bool,
    threads: Option<usize>,
    cancellation: &CancellationToken,
) -> Result<HashMap<String, TaxonId>> {
    let save_folder = save_folder.as_ref();
    fs::create_dir_all(save_folder)?;
    lookup_parsed_a2t(save_folder, accessions, low_memory, wgs, threads, cancellation)
}

/// Look up accessions assigned directly to the requested taxids without loading
/// the taxonomy tree.
pub fn lookup_taxid_accessions(
    save_folder: impl AsRef<Path>,
    taxa: &[TaxonId],
    low_memory: bool,
    wgs: bool,
    threads: Option<usize>,
) -> Result<HashSet<String>> {
    lookup_taxid_accessions_with_cancel(
        save_folder,
        taxa,
        low_memory,
        wgs,
        threads,
        &CancellationToken::default(),
    )
}

/// [`lookup_taxid_accessions`] that stops early when `cancellation` is triggered.
pub fn lookup_taxid_accessions_with_cancel(
    save_folder: impl AsRef<Path>,
    taxa: &[TaxonId],
    low_memory: bool,
    wgs: bool,
    threads: Option<usize>,
    cancellation: &CancellationToken,
) -> Result<HashSet<String>> {
    let save_folder = save_folder.as_ref();
    fs::create_dir_all(save_folder)?;
    lookup_t2a(save_folder, taxa, low_memory, wgs, threads, cancellation)
}

fn lookup_t2a(
    save_folder: &Path,
    taxa: &[TaxonId],
    low_memory: bool,
    wgs: bool,
    threads: Option<usize>,
    cancellation: &CancellationToken,
) -> Result<HashSet<String>> {
    let requested = taxa.iter().copied().collect::<HashSet<_>>();
    if requested.is_empty() {
        return Ok(HashSet::new());
    }
    if low_memory {
        let paths = ensure_a2t_files(save_folder, wgs)?;
        let decoder = lookup_decoder(crate::threads::resolve(threads)?)?;
        let partials = paths
            .par_iter()
            .map(|path| {
                let mut found = HashSet::new();
                scan_lookup_rows(path, &decoder, cancellation, |accession, taxid| {
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

    ensure_a2t_db(
        save_folder,
        false,
        wgs,
        crate::threads::resolve(threads)?,
        cancellation,
    )?;
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
    /// The installed file cannot be used as-is and must be built again.
    must_rebuild: bool,
}

/// Sizes and validators recorded for a source so an unchanged dump can be
/// recognised without downloading it again.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct SourceMeta {
    etag: Option<String>,
    last_modified: Option<String>,
    size: Option<i64>,
}

impl SourceMeta {
    /// Only treat a source as unchanged on positive evidence. A response that
    /// carried no validator at all matches nothing.
    fn matches(&self, other: &SourceMeta) -> bool {
        if self.etag.is_some() && self.etag == other.etag {
            return true;
        }
        self.size.is_some()
            && self.size == other.size
            && self.last_modified.is_some()
            && self.last_modified == other.last_modified
    }
}

fn header_value(response: &reqwest::blocking::Response, name: &str) -> Option<String> {
    response
        .headers()
        .get(name)
        .and_then(|value| value.to_str().ok())
        .map(str::to_owned)
}

fn response_meta(response: &reqwest::blocking::Response) -> SourceMeta {
    SourceMeta {
        etag: header_value(response, "etag"),
        last_modified: header_value(response, "last-modified"),
        size: header_value(response, "content-length").and_then(|value| value.parse().ok()),
    }
}

fn remote_source_meta(url: &str) -> Result<SourceMeta> {
    let response = Client::new().head(url).send()?.error_for_status()?;
    Ok(response_meta(&response))
}

/// Whether a source needs re-downloading, judged against what was recorded when
/// it was last fetched.
///
/// Being unable to reach the server is not evidence of a change: an offline run
/// keeps using the database it already has instead of starting a refresh it
/// cannot complete.
fn source_is_stale(connection: &Connection, source: &str) -> Result<bool> {
    let stored = stored_source_meta(connection, source)?;
    match remote_source_meta(&source_url(source)) {
        Err(_) => Ok(false),
        Ok(remote) => Ok(!stored.is_some_and(|stored| stored.matches(&remote))),
    }
}

fn stored_source_meta(connection: &Connection, source: &str) -> Result<Option<SourceMeta>> {
    let mut statement = connection
        .prepare("SELECT etag, last_modified, size FROM a2t_sources WHERE source = ?1")?;
    let mut rows = statement.query([source])?;
    let Some(row) = rows.next()? else {
        return Ok(None);
    };
    Ok(Some(SourceMeta {
        etag: row.get(0)?,
        last_modified: row.get(1)?,
        size: row.get(2)?,
    }))
}

fn record_source(
    connection: &Connection,
    source: &str,
    status: &str,
    meta: &SourceMeta,
    row_count: Option<u64>,
) -> Result<()> {
    connection.execute(
        "INSERT INTO a2t_sources(source, status, etag, last_modified, size, row_count)
         VALUES (?1, ?2, ?3, ?4, ?5, ?6)
         ON CONFLICT(source) DO UPDATE SET
             status = excluded.status,
             etag = excluded.etag,
             last_modified = excluded.last_modified,
             size = excluded.size,
             row_count = excluded.row_count",
        params![
            source,
            status,
            meta.etag,
            meta.last_modified,
            meta.size,
            row_count.map(|value| value as i64)
        ],
    )?;
    Ok(())
}

fn source_url(source: &str) -> String {
    format!("{A2T_BASE_URL}/{source}")
}

fn requested_a2t_sources(wgs: bool) -> Vec<&'static str> {
    let mut sources = vec![GB_FILE];
    if wgs {
        sources.push(WGS_FILE);
    }
    sources
}

/// Interrupts a connection's in-flight statement once `cancellation` trips.
///
/// `CREATE INDEX` over a billion rows is a single statement that runs for a long
/// time inside SQLite, so a cancellation flag checked between our own rows never
/// gets a turn. `InterruptHandle` is `Send`, so a watchdog thread can reach into
/// the running statement and make it return `SQLITE_INTERRUPT`.
struct SqlInterruptGuard {
    stop: Arc<AtomicBool>,
    handle: Option<thread::JoinHandle<()>>,
}

impl SqlInterruptGuard {
    fn new(connection: &Connection, cancellation: &CancellationToken) -> Self {
        let stop = Arc::new(AtomicBool::new(false));
        let interrupt = connection.get_interrupt_handle();
        let cancellation = cancellation.clone();
        let worker_stop = Arc::clone(&stop);
        let handle = thread::spawn(move || {
            while !worker_stop.load(AtomicOrdering::Relaxed) {
                if cancellation.is_cancelled() {
                    interrupt.interrupt();
                    return;
                }
                thread::sleep(Duration::from_millis(50));
            }
        });
        Self {
            stop,
            handle: Some(handle),
        }
    }
}

impl Drop for SqlInterruptGuard {
    fn drop(&mut self) {
        self.stop.store(true, AtomicOrdering::Relaxed);
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

type RowBatch = Vec<(Box<str>, TaxonId)>;

/// One decompressed source feeding the merge, read on its own thread.
struct MergeSource {
    receiver: Receiver<Result<RowBatch>>,
    handle: Option<thread::JoinHandle<()>>,
    batch: RowBatch,
    index: usize,
    finished: bool,
    rows: u64,
}

impl MergeSource {
    fn spawn(path: &Path, cancellation: &CancellationToken, threads: usize) -> Result<Self> {
        let (sender, receiver): (SyncSender<Result<RowBatch>>, _) = sync_channel(2);
        let path = path.to_path_buf();
        let cancellation = cancellation.clone();
        let handle = thread::spawn(move || {
            if let Err(error) = read_source_batches(&path, &cancellation, &sender, threads) {
                let _ = sender.send(Err(error));
            }
        });
        Ok(Self {
            receiver,
            handle: Some(handle),
            batch: Vec::new(),
            index: 0,
            finished: false,
            rows: 0,
        })
    }

    fn ensure_ready(&mut self) -> Result<()> {
        while !self.finished && self.index == self.batch.len() {
            match self.receiver.recv() {
                Ok(batch) => {
                    self.batch = batch?;
                    self.index = 0;
                }
                Err(_) => self.finished = true,
            }
        }
        Ok(())
    }

    fn current(&self) -> Option<&(Box<str>, TaxonId)> {
        self.batch.get(self.index)
    }

    fn take_current(&mut self) -> (Box<str>, TaxonId) {
        let row = std::mem::replace(&mut self.batch[self.index], (Box::from(""), 0));
        self.index += 1;
        self.rows += 1;
        row
    }
}

impl Drop for MergeSource {
    fn drop(&mut self) {
        // Dropping the receiver makes the reader's `send` fail, so it unwinds
        // rather than blocking forever on a full channel.
        if let Some(handle) = self.handle.take() {
            let _ = handle.join();
        }
    }
}

fn read_source_batches(
    path: &Path,
    cancellation: &CancellationToken,
    sender: &SyncSender<Result<RowBatch>>,
    threads: usize,
) -> Result<()> {
    let decoder = build_decoder(threads)?;
    let mut reader = BufReader::with_capacity(4 << 20, decoder.open(path)?);
    let mut batch: RowBatch = Vec::with_capacity(MERGE_BATCH_ROWS);
    let mut previous: Option<Box<str>> = None;
    let mut failed = false;
    let result = scan_reader(&mut reader, cancellation, &mut |accession, taxid| {
        // The merge and the incremental diff both assume each dump is already
        // sorted by accession in BINARY order, which NCBI's dumps are. Verify it
        // rather than silently producing a corrupt index if that ever changes.
        if let Some(previous) = &previous
            && accession < previous.as_ref()
        {
            bail!(
                "{} is not sorted by accession ({previous} precedes {accession})",
                path.display()
            );
        }
        previous = Some(Box::from(accession));
        batch.push((Box::from(accession), taxid));
        if batch.len() == MERGE_BATCH_ROWS {
            let full = std::mem::replace(&mut batch, Vec::with_capacity(MERGE_BATCH_ROWS));
            if sender.send(Ok(full)).is_err() {
                failed = true;
                return Ok(false);
            }
        }
        Ok(true)
    })
    .with_context(|| format!("failed to scan {}", path.display()));
    result?;
    if !failed && !batch.is_empty() {
        let _ = sender.send(Ok(batch));
    }
    Ok(())
}

/// Ascending merge of the per-source readers, with duplicate accessions
/// collapsed to their last occurrence.
struct MergeStream {
    sources: Vec<MergeSource>,
}

impl MergeStream {
    fn open(paths: &[PathBuf], cancellation: &CancellationToken, threads: usize) -> Result<Self> {
        // Every source is decompressed on its own reader thread at the same
        // time, so the budget is shared between them rather than applied to
        // each one.
        let per_source = (threads / paths.len().max(1)).max(1);
        let sources = paths
            .iter()
            .map(|path| MergeSource::spawn(path, cancellation, per_source))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self { sources })
    }

    fn lowest(&mut self) -> Result<Option<usize>> {
        for source in self.sources.iter_mut() {
            source.ensure_ready()?;
        }
        let mut best: Option<usize> = None;
        for (index, source) in self.sources.iter().enumerate() {
            let Some((accession, _)) = source.current() else {
                continue;
            };
            let better = match best {
                None => true,
                Some(current) => {
                    accession.as_ref() < self.sources[current].current().unwrap().0.as_ref()
                }
            };
            if better {
                best = Some(index);
            }
        }
        Ok(best)
    }

    fn next_row(&mut self) -> Result<Option<(Box<str>, TaxonId)>> {
        let Some(index) = self.lowest()? else {
            return Ok(None);
        };
        let mut row = self.sources[index].take_current();
        // An accession present in more than one dump, or repeated within one,
        // must not reach a UNIQUE key twice. Last occurrence wins.
        while let Some(next) = self.lowest()? {
            if self.sources[next].current().unwrap().0 != row.0 {
                break;
            }
            row = self.sources[next].take_current();
        }
        Ok(Some(row))
    }

    fn rows_per_source(&self) -> Vec<u64> {
        self.sources.iter().map(|source| source.rows).collect()
    }
}

fn schema_version(connection: &Connection) -> Result<Option<i64>> {
    let has_meta: bool = connection.query_row(
        "SELECT EXISTS(SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'a2t_meta')",
        [],
        |row| row.get(0),
    )?;
    if !has_meta {
        return Ok(None);
    }
    let mut statement =
        connection.prepare("SELECT value FROM a2t_meta WHERE key = 'schema_version'")?;
    let mut rows = statement.query([])?;
    let Some(row) = rows.next()? else {
        return Ok(None);
    };
    Ok(row.get::<_, String>(0)?.parse::<i64>().ok())
}

fn create_a2t_schema(connection: &Connection) -> Result<()> {
    connection.execute_batch(A2T_SCHEMA_SQL)?;
    connection.execute(
        "INSERT INTO a2t_meta(key, value) VALUES ('schema_version', ?1)",
        [SCHEMA_VERSION.to_string()],
    )?;
    Ok(())
}

fn inspect_a2t_database(db_path: &Path) -> Result<AccessionDatabaseState> {
    let rebuild = AccessionDatabaseState {
        loaded_sources: HashSet::new(),
        must_rebuild: true,
    };
    // Anything unreadable at the destination - truncated, corrupt, or not a
    // database at all - is replaced rather than reported as a failure.
    match read_a2t_database_state(db_path) {
        Ok(state) => Ok(state),
        Err(error) => {
            eprintln!(
                "taxutils: cannot read {} ({error:#}); rebuilding",
                db_path.display()
            );
            Ok(rebuild)
        }
    }
}

fn read_a2t_database_state(db_path: &Path) -> Result<AccessionDatabaseState> {
    let connection = Connection::open(db_path)?;
    let rebuild = AccessionDatabaseState {
        loaded_sources: HashSet::new(),
        must_rebuild: true,
    };

    match schema_version(&connection)? {
        Some(version) if version == SCHEMA_VERSION => {}
        Some(version) => {
            eprintln!(
                "taxutils: accession database is schema v{version}, this build needs \
                 v{SCHEMA_VERSION}; rebuilding once into the compact format"
            );
            return Ok(rebuild);
        }
        None => {
            eprintln!(
                "taxutils: accession database predates schema v{SCHEMA_VERSION}; \
                 rebuilding once into the compact format"
            );
            return Ok(rebuild);
        }
    }

    let incomplete: bool = connection.query_row(
        "SELECT EXISTS(SELECT 1 FROM a2t_sources WHERE status != 'complete')",
        [],
        |row| row.get(0),
    )?;
    let has_rows: bool =
        connection.query_row("SELECT EXISTS(SELECT 1 FROM a2t LIMIT 1)", [], |row| {
            row.get(0)
        })?;
    if incomplete || !has_rows {
        return Ok(rebuild);
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
        return Ok(rebuild);
    }
    Ok(AccessionDatabaseState {
        loaded_sources,
        must_rebuild: false,
    })
}

fn configure_bulk_load(connection: &Connection, threads: usize) -> Result<()> {
    // Larger pages cut per-page overhead across a table of this size. It must be
    // set before anything is written, so this runs before the schema is created.
    connection.pragma_update(None, "page_size", 8192_i64)?;
    connection.pragma_update(None, "journal_mode", "OFF")?;
    connection.pragma_update(None, "synchronous", "OFF")?;
    connection.pragma_update(None, "locking_mode", "EXCLUSIVE")?;
    connection.pragma_update(None, "cache_size", -524_288_i64)?;
    connection.pragma_update(None, "temp_store", 1_i64)?;
    connection.pragma_update(None, "threads", threads as i64)?;
    Ok(())
}

/// Durable settings for mutating a database that is already installed.
fn configure_incremental(connection: &Connection, threads: usize) -> Result<()> {
    connection.pragma_update(None, "journal_mode", "WAL")?;
    connection.pragma_update(None, "synchronous", "NORMAL")?;
    connection.pragma_update(None, "cache_size", -262_144_i64)?;
    connection.pragma_update(None, "temp_store", 1_i64)?;
    connection.pragma_update(None, "threads", threads as i64)?;
    Ok(())
}

fn finish_bulk_load(connection: &Connection, cancellation: &CancellationToken) -> Result<()> {
    {
        let _guard = SqlInterruptGuard::new(connection, cancellation);
        // The only remaining sort. The table itself arrived in key order, so
        // there is no separate accession index to build.
        connection.execute_batch("CREATE INDEX idx_taxid ON a2t(taxid);")?;
    }
    cancellation.check_cancelled()?;
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

/// Load every source into a fresh database in one merged, ascending pass.
fn build_a2t_database_atomic(
    sources: &[PathBuf],
    metas: &HashMap<String, SourceMeta>,
    db_path: &Path,
    cancellation: &CancellationToken,
    threads: usize,
) -> Result<()> {
    let parent = db_path.parent().unwrap_or_else(|| Path::new("."));
    let temporary = tempfile::NamedTempFile::new_in(parent)?.into_temp_path();
    let mut connection = Connection::open(&temporary)?;
    configure_bulk_load(&connection, threads)?;
    create_a2t_schema(&connection)?;

    for source in sources {
        let name = source_name(source);
        record_source(&connection, &name, "loading", &SourceMeta::default(), None)?;
    }

    let rows_per_source = {
        let mut merge = MergeStream::open(sources, cancellation, threads)?;
        let transaction = connection.transaction()?;
        {
            // Keys arrive in primary-key order, so each page fills once and is
            // never revisited. OR REPLACE covers an accession seen in two dumps.
            let mut insert = transaction
                .prepare("INSERT OR REPLACE INTO a2t(accession, taxid) VALUES (?1, ?2)")?;
            while let Some((accession, taxid)) = merge.next_row()? {
                insert.execute(params![accession.as_ref(), taxid])?;
            }
        }
        transaction.commit()?;
        merge.rows_per_source()
    };
    cancellation.check_cancelled()?;

    for (source, rows) in sources.iter().zip(rows_per_source) {
        let name = source_name(source);
        let meta = metas.get(&name).cloned().unwrap_or_default();
        record_source(&connection, &name, "complete", &meta, Some(rows))?;
    }

    finish_bulk_load(&connection, cancellation)?;
    connection.close().map_err(|(_, error)| error)?;
    install_temporary_database(temporary, db_path)
}

#[derive(Debug, Default)]
struct RefreshStats {
    inserted: u64,
    updated: u64,
    deleted: u64,
}

/// Bring an installed database up to date without rebuilding it.
///
/// Both the merged input and `a2t` are ordered by accession, so one lockstep
/// pass classifies every row as an insert, an update, or a deletion. The diff is
/// staged in a temp table because mutating `a2t` while a cursor walks it is not
/// reliable.
fn refresh_a2t_database(
    db_path: &Path,
    sources: &[PathBuf],
    metas: &HashMap<String, SourceMeta>,
    cancellation: &CancellationToken,
    threads: usize,
) -> Result<RefreshStats> {
    let mut connection = Connection::open(db_path)?;
    configure_incremental(&connection, threads)?;
    connection.execute_batch(
        "DROP TABLE IF EXISTS temp.a2t_delta;
         CREATE TEMP TABLE a2t_delta (accession TEXT PRIMARY KEY, taxid INTEGER, op INTEGER)
             WITHOUT ROWID;",
    )?;

    let mut stats = RefreshStats::default();
    let rows_per_source = {
        let mut merge = MergeStream::open(sources, cancellation, threads)?;
        let mut table =
            connection.prepare("SELECT accession, taxid FROM a2t ORDER BY accession")?;
        let mut rows = table.query([])?;
        let mut existing: Option<(String, TaxonId)> = match rows.next()? {
            Some(row) => Some((row.get(0)?, row.get(1)?)),
            None => None,
        };
        let mut incoming = merge.next_row()?;

        let delta = connection.unchecked_transaction()?;
        {
            let mut stage = delta
                .prepare("INSERT INTO temp.a2t_delta(accession, taxid, op) VALUES (?1, ?2, ?3)")?;
            let mut scanned: u64 = 0;
            loop {
                scanned += 1;
                if scanned.is_multiple_of(CANCEL_CHECK_ROWS) {
                    cancellation.check_cancelled()?;
                }
                match (&incoming, &existing) {
                    (None, None) => break,
                    (Some((accession, taxid)), None) => {
                        stage.execute(params![accession.as_ref(), taxid, OP_INSERT])?;
                        stats.inserted += 1;
                        incoming = merge.next_row()?;
                    }
                    (None, Some((accession, _))) => {
                        stage.execute(params![accession, None::<TaxonId>, OP_DELETE])?;
                        stats.deleted += 1;
                        existing = match rows.next()? {
                            Some(row) => Some((row.get(0)?, row.get(1)?)),
                            None => None,
                        };
                    }
                    (Some((new_accession, new_taxid)), Some((old_accession, old_taxid))) => {
                        match new_accession.as_ref().cmp(old_accession.as_str()) {
                            Ordering::Less => {
                                stage.execute(params![
                                    new_accession.as_ref(),
                                    new_taxid,
                                    OP_INSERT
                                ])?;
                                stats.inserted += 1;
                                incoming = merge.next_row()?;
                            }
                            Ordering::Greater => {
                                stage.execute(params![
                                    old_accession,
                                    None::<TaxonId>,
                                    OP_DELETE
                                ])?;
                                stats.deleted += 1;
                                existing = match rows.next()? {
                                    Some(row) => Some((row.get(0)?, row.get(1)?)),
                                    None => None,
                                };
                            }
                            Ordering::Equal => {
                                if new_taxid != old_taxid {
                                    stage.execute(params![
                                        new_accession.as_ref(),
                                        new_taxid,
                                        OP_UPDATE
                                    ])?;
                                    stats.updated += 1;
                                }
                                incoming = merge.next_row()?;
                                existing = match rows.next()? {
                                    Some(row) => Some((row.get(0)?, row.get(1)?)),
                                    None => None,
                                };
                            }
                        }
                    }
                }
            }
        }
        delta.commit()?;
        merge.rows_per_source()
    };
    cancellation.check_cancelled()?;

    // Apply and record the new source validators together, so a crash cannot
    // leave the metadata claiming a refresh that did not land.
    let transaction = connection.transaction()?;
    {
        let _guard = SqlInterruptGuard::new(&transaction, cancellation);
        transaction.execute_batch(
            "INSERT OR REPLACE INTO a2t(accession, taxid)
                 SELECT accession, taxid FROM temp.a2t_delta WHERE op IN (0, 1);
             DELETE FROM a2t WHERE accession IN
                 (SELECT accession FROM temp.a2t_delta WHERE op = 2);",
        )?;
    }
    for (source, rows) in sources.iter().zip(rows_per_source) {
        let name = source_name(source);
        let meta = metas.get(&name).cloned().unwrap_or_default();
        record_source(&transaction, &name, "complete", &meta, Some(rows))?;
    }
    transaction.commit()?;
    connection.execute_batch("DROP TABLE IF EXISTS temp.a2t_delta;")?;
    connection.pragma_update(None, "journal_mode", "DELETE")?;
    connection.pragma_update(None, "synchronous", "FULL")?;
    Ok(stats)
}

const OP_INSERT: i64 = 0;
const OP_UPDATE: i64 = 1;
const OP_DELETE: i64 = 2;

fn source_name(path: &Path) -> String {
    path.file_name().unwrap().to_string_lossy().into_owned()
}

/// Remove temporary databases orphaned by a hard kill.
///
/// A cancelled build cleans up after itself through `TempPath`, but `SIGKILL`
/// leaves multi-gigabyte files behind in the save folder.
fn sweep_stale_temporaries(save_folder: &Path) -> Result<()> {
    let Ok(entries) = fs::read_dir(save_folder) else {
        return Ok(());
    };
    for entry in entries.flatten() {
        let name = entry.file_name();
        let Some(name) = name.to_str() else { continue };
        if !name.starts_with(".tmp") {
            continue;
        }
        if entry.path().is_file() {
            let _ = fs::remove_file(entry.path());
        }
    }
    Ok(())
}

/// Collect the current validators for the sources that will be merged, and
/// download whatever is missing or has changed upstream.
fn prepare_sources(
    save_folder: &Path,
    names: &[String],
    stale: &HashSet<String>,
    previous: &HashMap<String, SourceMeta>,
) -> Result<(Vec<PathBuf>, HashMap<String, SourceMeta>)> {
    let paths = names
        .iter()
        .map(|name| save_folder.join(name))
        .collect::<Vec<_>>();
    let metas = names
        .par_iter()
        .map(|name| {
            let path = save_folder.join(name);
            if stale.contains(name) || !path.exists() {
                let meta = download_file_with_meta(&source_url(name), &path)?;
                return Ok((name.clone(), meta));
            }
            Ok((
                name.clone(),
                previous.get(name).cloned().unwrap_or_default(),
            ))
        })
        .collect::<Result<Vec<_>>>()?
        .into_iter()
        .collect();
    Ok((paths, metas))
}

fn ensure_a2t_db(
    save_folder: &Path,
    refresh: bool,
    wgs: bool,
    threads: usize,
    cancellation: &CancellationToken,
) -> Result<()> {
    fs::create_dir_all(save_folder)?;
    sweep_stale_temporaries(save_folder)?;
    let db_path = save_folder.join(DB_FILE);
    let requested = requested_a2t_sources(wgs)
        .into_iter()
        .map(str::to_owned)
        .collect::<Vec<_>>();

    if db_path.exists() {
        let state = inspect_a2t_database(&db_path)?;
        if !state.must_rebuild {
            let connection = Connection::open(&db_path)?;
            let missing = requested
                .iter()
                .filter(|source| !state.loaded_sources.contains(*source))
                .cloned()
                .collect::<HashSet<_>>();

            // Anything already loaded stays in the merge set. Otherwise a run
            // with `wgs = false` against a database built with WGS would read
            // every WGS accession as withdrawn and delete it.
            let mut merge_names = state.loaded_sources.iter().cloned().collect::<HashSet<_>>();
            merge_names.extend(requested.iter().cloned());
            let mut merge_names = merge_names.into_iter().collect::<Vec<_>>();
            merge_names.sort();

            let mut stale = missing.clone();
            let mut previous = HashMap::new();
            for source in &merge_names {
                if let Some(meta) = stored_source_meta(&connection, source)? {
                    previous.insert(source.clone(), meta);
                }
                // Only a caller that asked for a refresh pays for a round trip to
                // NCBI. Simply opening the database stays offline.
                if refresh && !stale.contains(source) && source_is_stale(&connection, source)? {
                    stale.insert(source.clone());
                }
            }
            drop(connection);

            if stale.is_empty() {
                return Ok(());
            }

            let (paths, metas) = prepare_sources(save_folder, &merge_names, &stale, &previous)?;
            let stats =
                refresh_a2t_database(&db_path, &paths, &metas, cancellation, threads)?;
            eprintln!(
                "taxutils: refreshed accession database ({} inserted, {} updated, {} deleted)",
                stats.inserted, stats.updated, stats.deleted
            );
            return Ok(());
        }
    }

    // Nothing usable is installed: download whatever is missing and build.
    let (paths, metas) =
        prepare_sources(save_folder, &requested, &HashSet::new(), &HashMap::new())?;
    build_a2t_database_atomic(&paths, &metas, &db_path, cancellation, threads)?;
    Ok(())
}

/// Build, refresh, or validate the shared SQLite accession index.
///
/// A new database is assembled beside the destination and atomically installed
/// only once every row and index is complete. An existing database is brought up
/// to date in place, row by row, rather than rebuilt.
pub fn ensure_accession_database(
    save_folder: impl AsRef<Path>,
    options: AccessionDatabaseOptions,
) -> Result<PathBuf> {
    ensure_accession_database_with_cancel(save_folder, options, &CancellationToken::default())
}

/// [`ensure_accession_database`] that stops early when `cancellation` is triggered.
///
/// A cancelled build leaves no temporary file behind and leaves any database
/// already installed at the destination untouched.
pub fn ensure_accession_database_with_cancel(
    save_folder: impl AsRef<Path>,
    options: AccessionDatabaseOptions,
    cancellation: &CancellationToken,
) -> Result<PathBuf> {
    let save_folder = save_folder.as_ref();
    ensure_a2t_db(
        save_folder,
        options.refresh,
        options.wgs,
        crate::threads::resolve(options.threads)?,
        cancellation,
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
            scan_reader(
                &mut reader,
                &CancellationToken::default(),
                &mut |accession: &str, taxid| {
                    rows.push((accession.to_owned(), taxid));
                    Ok(true)
                },
            )
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
                lookup_accession_taxids(dir.path(), HashSet::new(), low_memory, true, None)
                    .unwrap()
                    .is_empty()
            );
            assert!(
                lookup_taxid_accessions(dir.path(), &[], low_memory, true, None)
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
                        lookup_accession_taxids(dir.path(), queries.clone(), true, false, None).unwrap();
                    assert_eq!(
                        gb,
                        HashMap::from([
                            ("NC_000001.1".to_owned(), 13),
                            ("NC_000002.1".to_owned(), 15)
                        ])
                    );
                    let both =
                        lookup_accession_taxids(dir.path(), queries.clone(), true, true, None).unwrap();
                    assert_eq!(both["NC_000001.1"], 99); // WGS keeps its existing precedence.
                    assert_eq!(both["ABCD01000001.1"], 15);
                    assert_eq!(
                        lookup_taxid_accessions(dir.path(), &[15, 15], true, false, None).unwrap(),
                        HashSet::from(["NC_000002.1".to_owned()])
                    );
                    assert_eq!(
                        lookup_taxid_accessions(dir.path(), &[15], true, true, None).unwrap(),
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
                , None)
                .is_err()
            );
            assert!(lookup_taxid_accessions(dir.path(), &[13], true, false, None).is_err());
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
            None,
        )
            .unwrap();
            assert_eq!(found["NC_000000000.1"], 0);
            assert_eq!(found["NC_000099999.1"], 999);
            assert!(lookup_taxid_accessions(dir.path(), &[999], true, false, None).is_err());
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
        let token = CancellationToken::default();
        let scanned = lookup_a2t(dir.path(), &accessions, true, true, None, &token).unwrap();
        let indexed = lookup_a2t(dir.path(), &accessions, false, true, None, &token).unwrap();
        let direct = lookup_accession_taxids(
            dir.path(),
            accessions.iter().cloned().collect(),
            false,
            true,
        None,
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
        // `a2t` is keyed on the accession itself, so the forward lookup seeks the
        // table's own primary-key index and there is no second copy to maintain.
        assert!(
            details
                .iter()
                .any(|detail| detail.contains("SEARCH a USING PRIMARY KEY")),
            "forward lookup did not seek the accession primary key: {details:?}"
        );
        assert!(
            !details
                .iter()
                .any(|detail| detail == "SCAN a" || detail.starts_with("SCAN a ")),
            "forward lookup scanned the accession table: {details:?}"
        );

        let scanned_reverse = lookup_t2a(dir.path(), &[20, 30], true, true, None, &token).unwrap();
        let indexed_reverse = lookup_t2a(dir.path(), &[20, 30], false, true, None, &token).unwrap();
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
        // On a WITHOUT ROWID table the index carries the primary key, so
        // idx_taxid is covering: the reverse lookup is answered from the index
        // alone and never visits the table.
        assert!(
            details.iter().any(|detail| {
                detail.contains("SEARCH a USING COVERING INDEX idx_taxid")
                    && detail.contains("taxid=?")
            }),
            "reverse lookup must search idx_taxid as a covering index: {details:?}"
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
            let scanned = lookup_taxid_accessions(dir.path(), &taxa, true, false, None).unwrap();
            let indexed = lookup_taxid_accessions(dir.path(), &taxa, false, false, None).unwrap();
            assert_eq!(indexed, scanned, "reverse lookup mismatch for {taxa:?}");
        }
        assert_eq!(
            lookup_taxid_accessions(dir.path(), &[12059], false, false, None).unwrap(),
            HashSet::from(["NC_000001.1".to_owned(), "NC_000002.1".to_owned()])
        );
    }

    #[test]
    fn database_build_is_atomic_and_retains_sources() {
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join(GB_FILE);
        write_a2t_fixture(&source, &[("NC_000001.1", 10), ("NC_000002.1", 20)]);

        let db_path =
            ensure_accession_database(dir.path(), AccessionDatabaseOptions::default()).unwrap();
        assert!(db_path.exists());
        // The compressed source is kept: low-memory lookups read it directly.
        assert!(source.exists());
        assert!(
            fs::read_dir(dir.path()).unwrap().all(|entry| !entry
                .unwrap()
                .file_name()
                .to_string_lossy()
                .contains("tmp"))
        );

        // A complete database is reused as-is on the next call.
        ensure_accession_database(dir.path(), AccessionDatabaseOptions::default()).unwrap();
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
        // The accession index is gone: the WITHOUT ROWID table serves that role.
        assert_eq!(indexes, HashSet::from(["idx_taxid".to_owned()]));
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

        assert!(
            build_a2t_database_atomic(
                &[invalid_source],
                &HashMap::new(),
                &db_path,
                &CancellationToken::default(),
            1,
        )
            .is_err()
        );
        assert_eq!(fs::read(&db_path).unwrap(), b"existing database sentinel");
    }

    #[test]
    fn empty_database_with_complete_metadata_is_rebuilt() {
        let dir = tempfile::tempdir().unwrap();
        let db_path = dir.path().join(DB_FILE);
        let connection = Connection::open(&db_path).unwrap();
        create_a2t_schema(&connection).unwrap();
        connection
            .execute(
                "INSERT INTO a2t_sources(source, status) VALUES (?1, 'complete')",
                [GB_FILE],
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

    /// A database written by an older release has no schema marker, so it must
    /// be rebuilt rather than read with the wrong layout.
    #[test]
    fn legacy_schema_database_is_rebuilt() {
        let dir = tempfile::tempdir().unwrap();
        let db_path = dir.path().join(DB_FILE);
        let connection = Connection::open(&db_path).unwrap();
        connection
            .execute_batch(
                "CREATE TABLE a2t (accession TEXT, taxid INTEGER);
                 CREATE TABLE a2t_sources (source TEXT PRIMARY KEY, status TEXT);
                 INSERT INTO a2t VALUES ('NC_000001.1', 999);
                 INSERT INTO a2t_sources VALUES ('nucl_gb.accession2taxid.gz', 'complete');
                 CREATE INDEX idx_accession ON a2t(accession);
                 CREATE INDEX idx_taxid ON a2t(taxid);",
            )
            .unwrap();
        drop(connection);
        assert!(inspect_a2t_database(&db_path).unwrap().must_rebuild);

        write_a2t_fixture(&dir.path().join(GB_FILE), &[("NC_000001.1", 10)]);
        ensure_accession_database(dir.path(), AccessionDatabaseOptions::default()).unwrap();

        let connection = Connection::open(&db_path).unwrap();
        assert_eq!(schema_version(&connection).unwrap(), Some(SCHEMA_VERSION));
        let taxid: TaxonId = connection
            .query_row(
                "SELECT taxid FROM a2t WHERE accession = 'NC_000001.1'",
                [],
                |row| row.get(0),
            )
            .unwrap();
        assert_eq!(taxid, 10, "stale row from the legacy database survived");
        let indexes = connection
            .prepare("SELECT name FROM sqlite_master WHERE type = 'index' AND name LIKE 'idx_%'")
            .unwrap()
            .query_map([], |row| row.get::<_, String>(0))
            .unwrap()
            .collect::<rusqlite::Result<HashSet<_>>>()
            .unwrap();
        assert_eq!(indexes, HashSet::from(["idx_taxid".to_owned()]));
    }

    /// Sources are merged into one ascending stream, and the two NCBI dumps
    /// interleave rather than sitting in disjoint key ranges.
    #[test]
    fn merge_interleaves_sources_in_key_order_and_collapses_duplicates() {
        let dir = tempfile::tempdir().unwrap();
        write_a2t_fixture(
            &dir.path().join(GB_FILE),
            &[("A00001.1", 1), ("AB000001.1", 2), ("NC_000001.1", 3)],
        );
        write_a2t_fixture(
            &dir.path().join(WGS_FILE),
            &[
                ("AAAA01000001.1", 4),
                ("AB000001.1", 5),
                ("ZZZZ01000001.1", 6),
            ],
        );
        let paths = vec![dir.path().join(GB_FILE), dir.path().join(WGS_FILE)];
        let token = CancellationToken::default();
        let mut merge = MergeStream::open(&paths, &token, 1).unwrap();
        let mut rows = Vec::new();
        while let Some((accession, taxid)) = merge.next_row().unwrap() {
            rows.push((accession.to_string(), taxid));
        }
        assert_eq!(
            rows,
            vec![
                ("A00001.1".to_owned(), 1),
                ("AAAA01000001.1".to_owned(), 4),
                // Present in both dumps; the later source wins, and it is emitted once.
                ("AB000001.1".to_owned(), 5),
                ("NC_000001.1".to_owned(), 3),
                ("ZZZZ01000001.1".to_owned(), 6),
            ]
        );
    }

    /// An out-of-order dump would silently corrupt both the merge and the diff,
    /// so it must be rejected instead of trusted.
    #[test]
    fn unsorted_source_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join(GB_FILE);
        write_a2t_fixture(&source, &[("NC_000002.1", 10), ("NC_000001.1", 20)]);
        let token = CancellationToken::default();
        let mut merge = MergeStream::open(&[source], &token, 1).unwrap();
        let error = loop {
            match merge.next_row() {
                Ok(Some(_)) => continue,
                Ok(None) => panic!("unsorted source was accepted"),
                Err(error) => break error,
            }
        };
        assert!(
            format!("{error:#}").contains("not sorted"),
            "unexpected error: {error:#}"
        );
    }

    /// The refresh path must insert, update and delete without rebuilding.
    #[test]
    fn refresh_applies_row_level_difference() {
        let dir = tempfile::tempdir().unwrap();
        let source = dir.path().join(GB_FILE);
        write_a2t_fixture(
            &source,
            &[
                ("NC_000001.1", 10),
                ("NC_000002.1", 20),
                ("NC_000003.1", 30),
            ],
        );
        let db_path =
            ensure_accession_database(dir.path(), AccessionDatabaseOptions::default()).unwrap();
        let before = fs::metadata(&db_path).unwrap().modified().unwrap();

        // NC_000002 changes taxid, NC_000003 is withdrawn, NC_000004 is new.
        write_a2t_fixture(
            &source,
            &[
                ("NC_000001.1", 10),
                ("NC_000002.1", 99),
                ("NC_000004.1", 40),
            ],
        );
        let metas = HashMap::new();
        let stats =
            refresh_a2t_database(&db_path, &[source], &metas, &CancellationToken::default(), 1)
                .unwrap();
        assert_eq!(stats.inserted, 1);
        assert_eq!(stats.updated, 1);
        assert_eq!(stats.deleted, 1);

        let connection = Connection::open(&db_path).unwrap();
        let mut statement = connection
            .prepare("SELECT accession, taxid FROM a2t ORDER BY accession")
            .unwrap();
        let rows = statement
            .query_map([], |row| {
                Ok((row.get::<_, String>(0)?, row.get::<_, TaxonId>(1)?))
            })
            .unwrap()
            .collect::<rusqlite::Result<Vec<_>>>()
            .unwrap();
        assert_eq!(
            rows,
            vec![
                ("NC_000001.1".to_owned(), 10),
                ("NC_000002.1".to_owned(), 99),
                ("NC_000004.1".to_owned(), 40),
            ]
        );
        let _ = before;
    }

    /// Cancelling a build must leave the destination and the save folder as they
    /// were, with no partial database left lying around.
    #[test]
    fn cancelled_build_leaves_no_temporary_and_preserves_destination() {
        let dir = tempfile::tempdir().unwrap();
        let db_path = dir.path().join(DB_FILE);
        fs::write(&db_path, b"existing database sentinel").unwrap();
        let source = dir.path().join(GB_FILE);
        let rows = (0..200_000)
            .map(|index| (format!("NC_{index:09}.1"), 10))
            .collect::<Vec<_>>();
        write_a2t_fixture(
            &source,
            &rows
                .iter()
                .map(|(accession, taxid)| (accession.as_str(), *taxid))
                .collect::<Vec<_>>(),
        );

        let cancellation = CancellationToken::default();
        cancellation.cancel();
        let result = build_a2t_database_atomic(&[source], &HashMap::new(), &db_path, &cancellation, 1);
        assert!(result.is_err(), "cancelled build reported success");
        assert_eq!(fs::read(&db_path).unwrap(), b"existing database sentinel");
        let leftovers = fs::read_dir(dir.path())
            .unwrap()
            .filter_map(|entry| entry.ok())
            .filter(|entry| entry.file_name().to_string_lossy().starts_with(".tmp"))
            .count();
        assert_eq!(leftovers, 0, "cancelled build left a temporary database");
    }
}
