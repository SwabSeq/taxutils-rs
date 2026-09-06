use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use crate::TaxutilsOptions;
use crate::accession::{parse_accession, parse_accessions};
use crate::resources::AccessionTaxidIndex;

const IO_BUFFER_BYTES: usize = 1 << 20;
const CLEAN_BATCH_BYTES: usize = 8 << 20;

/// Cooperative cancellation shared with language bindings.
#[derive(Clone, Debug, Default)]
pub struct CancellationToken(Arc<AtomicBool>);

impl CancellationToken {
    pub fn cancel(&self) {
        self.0.store(true, Ordering::Relaxed);
    }

    pub fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }

    fn check(&self) -> Result<()> {
        self.check_cancelled()
    }

    /// Error out if cancellation has been requested.
    pub fn check_cancelled(&self) -> Result<()> {
        if self.is_cancelled() {
            bail!("operation cancelled");
        }
        Ok(())
    }
}

#[derive(Debug)]
struct FastaRecord {
    data: Vec<u8>,
    header_len: usize,
}

impl FastaRecord {
    fn header(&self) -> &[u8] {
        &self.data[..self.header_len]
    }

    fn write_to(&self, writer: &mut impl Write) -> Result<()> {
        writer.write_all(&self.data)?;
        Ok(())
    }
}

struct FastaReader<R> {
    reader: R,
    pending_header: Option<Vec<u8>>,
    finished: bool,
}

impl<R: BufRead> FastaReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            pending_header: None,
            finished: false,
        }
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = std::io::Result<FastaRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        let header = if let Some(header) = self.pending_header.take() {
            header
        } else {
            loop {
                let mut line = Vec::new();
                match self.reader.read_until(b'\n', &mut line) {
                    Ok(0) => {
                        self.finished = true;
                        return None;
                    }
                    Ok(_) if line.starts_with(b">") => break line,
                    Ok(_) => {}
                    Err(error) => return Some(Err(error)),
                }
            }
        };

        let header_len = header.len();
        let mut data = header;
        loop {
            let line_start = data.len();
            match self.reader.read_until(b'\n', &mut data) {
                Ok(0) => {
                    self.finished = true;
                    break;
                }
                Ok(_) if data[line_start] == b'>' => {
                    self.pending_header = Some(data.split_off(line_start));
                    break;
                }
                Ok(_) => {}
                Err(error) => return Some(Err(error)),
            }
        }
        Some(Ok(FastaRecord { data, header_len }))
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FilterMode {
    Keep,
    Remove,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FilterStats {
    pub kept: usize,
    pub removed: usize,
    pub missing_accession: usize,
    pub missing_taxid: usize,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GrepStats {
    pub requested: usize,
    pub scanned: usize,
    pub matched: usize,
    pub missing_accession: usize,
}

fn record_batch<R: BufRead>(
    records: &mut FastaReader<R>,
    max_records: usize,
    max_bytes: usize,
) -> Result<Vec<FastaRecord>> {
    let mut batch = Vec::new();
    let mut bytes = 0;
    while batch.len() < max_records && (bytes < max_bytes || batch.is_empty()) {
        let Some(record) = records.next() else {
            break;
        };
        let record = record?;
        bytes += record.data.len();
        batch.push(record);
    }
    Ok(batch)
}

fn write_accession_batch(
    output: &mut impl Write,
    headers: &mut Vec<(usize, Vec<u8>)>,
) -> Result<usize> {
    let accessions = headers
        .par_iter()
        .map(|(_, header)| parse_accession(&String::from_utf8_lossy(header), true))
        .collect::<Vec<_>>();
    for ((line_number, header), accession) in headers.iter().zip(&accessions) {
        if accession == "NA" {
            bail!(
                "No accession found in FASTA header on line {line_number}: {}",
                String::from_utf8_lossy(header).trim()
            );
        }
        writeln!(output, "{accession}")?;
    }
    let count = headers.len();
    headers.clear();
    Ok(count)
}

pub fn extract_accessions(
    fasta_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    batch_size: usize,
) -> Result<usize> {
    extract_accessions_with_cancel(
        fasta_path,
        output_path,
        batch_size,
        &CancellationToken::default(),
    )
}

pub fn extract_accessions_with_cancel(
    fasta_path: impl AsRef<Path>,
    output_path: impl AsRef<Path>,
    batch_size: usize,
    cancellation: &CancellationToken,
) -> Result<usize> {
    if batch_size < 1 {
        bail!("--batch-size must be at least 1");
    }
    let output_path = output_path.as_ref();
    let output_dir = output_path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_dir)?;
    let temporary = tempfile::NamedTempFile::new_in(output_dir)?;
    let mut output = BufWriter::with_capacity(IO_BUFFER_BYTES, temporary);
    let mut reader = BufReader::with_capacity(IO_BUFFER_BYTES, File::open(fasta_path)?);
    let mut headers = Vec::with_capacity(batch_size);
    let mut line = Vec::new();
    let mut count = 0;
    let mut line_number = 0;
    loop {
        cancellation.check()?;
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        line_number += 1;
        if !line.starts_with(b">") {
            continue;
        }
        headers.push((line_number, std::mem::take(&mut line)));
        if headers.len() == batch_size {
            count += write_accession_batch(&mut output, &mut headers)?;
        }
    }
    count += write_accession_batch(&mut output, &mut headers)?;
    output.flush()?;
    let temporary = output.into_inner().map_err(|error| error.into_error())?;
    temporary
        .persist(output_path)
        .map_err(|error| error.error)?;
    Ok(count)
}

#[derive(Clone, Copy, Debug)]
struct CleanHeader {
    line_number: usize,
    start: usize,
    end: usize,
}

fn write_clean_batch(
    output: &mut impl Write,
    data: &mut Vec<u8>,
    headers: &mut Vec<CleanHeader>,
    verbose: bool,
) -> Result<()> {
    let accessions = headers
        .par_iter()
        .map(|header| {
            parse_accession(
                &String::from_utf8_lossy(&data[header.start..header.end]),
                true,
            )
        })
        .collect::<Vec<_>>();
    let mut cursor = 0;
    for (header, accession) in headers.iter().zip(accessions) {
        output.write_all(&data[cursor..header.start])?;
        let original = &data[header.start..header.end];
        if accession == "NA" {
            if verbose {
                println!(
                    "NA accession line {}: {}",
                    header.line_number,
                    String::from_utf8_lossy(original).trim()
                );
            }
            bail!(
                "No accession found in FASTA header on line {}: {}",
                header.line_number,
                String::from_utf8_lossy(original).trim()
            );
        }
        writeln!(output, ">{accession}")?;
        cursor = header.end;
    }
    output.write_all(&data[cursor..])?;
    data.clear();
    headers.clear();
    Ok(())
}

pub fn clean_fasta_headers(
    input_path: impl AsRef<Path>,
    output_path: Option<&Path>,
    verbose: bool,
) -> Result<()> {
    clean_fasta_headers_with_cancel(
        input_path,
        output_path,
        verbose,
        &CancellationToken::default(),
    )
}

pub fn clean_fasta_headers_with_cancel(
    input_path: impl AsRef<Path>,
    output_path: Option<&Path>,
    verbose: bool,
    cancellation: &CancellationToken,
) -> Result<()> {
    let input_path = input_path.as_ref();
    let destination = output_path.unwrap_or(input_path);
    let output_dir = destination.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_dir)?;
    let temporary = tempfile::NamedTempFile::new_in(output_dir)?;
    let mut output = BufWriter::with_capacity(IO_BUFFER_BYTES, temporary);
    let mut reader = BufReader::with_capacity(IO_BUFFER_BYTES, File::open(input_path)?);
    let mut data = Vec::with_capacity(CLEAN_BATCH_BYTES);
    let mut headers = Vec::new();
    let mut line_number = 0;
    loop {
        cancellation.check()?;
        let mut reached_eof = false;
        while data.len() < CLEAN_BATCH_BYTES {
            let start = data.len();
            if reader.read_until(b'\n', &mut data)? == 0 {
                reached_eof = true;
                break;
            }
            line_number += 1;
            if data[start] == b'>' {
                headers.push(CleanHeader {
                    line_number,
                    start,
                    end: data.len(),
                });
            }
        }
        if data.is_empty() {
            break;
        }
        write_clean_batch(&mut output, &mut data, &mut headers, verbose)?;
        if reached_eof {
            break;
        }
    }
    output.flush()?;
    let temporary = output.into_inner().map_err(|error| error.into_error())?;
    temporary
        .persist(destination)
        .map_err(|error| error.error)?;
    Ok(())
}

fn read_query(value: &str) -> Result<String> {
    if Path::new(value).exists() {
        Ok(fs::read_to_string(value)?)
    } else {
        Ok(value.to_owned())
    }
}

fn grep_batch(
    output: &mut impl Write,
    batch: &[FastaRecord],
    requested: &HashSet<String>,
    version: bool,
    verbose: bool,
    stats: &mut GrepStats,
) -> Result<()> {
    let accessions = batch
        .par_iter()
        .map(|record| parse_accession(&String::from_utf8_lossy(record.header()), version))
        .collect::<Vec<_>>();
    for (record, accession) in batch.iter().zip(accessions) {
        stats.scanned += 1;
        if accession == "NA" {
            stats.missing_accession += 1;
            if verbose {
                println!(
                    "NA accession: {}",
                    String::from_utf8_lossy(record.header()).trim()
                );
            }
        } else if requested.contains(&accession) {
            record.write_to(output)?;
            stats.matched += 1;
        }
    }
    Ok(())
}

pub fn grep_fasta(
    input_path: impl AsRef<Path>,
    accession_query: &str,
    output_path: impl AsRef<Path>,
    version: bool,
    batch_size: usize,
    verbose: bool,
) -> Result<GrepStats> {
    grep_fasta_with_cancel(
        input_path,
        accession_query,
        output_path,
        version,
        batch_size,
        verbose,
        &CancellationToken::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn grep_fasta_with_cancel(
    input_path: impl AsRef<Path>,
    accession_query: &str,
    output_path: impl AsRef<Path>,
    version: bool,
    batch_size: usize,
    verbose: bool,
    cancellation: &CancellationToken,
) -> Result<GrepStats> {
    if batch_size < 1 {
        bail!("--batch-size must be at least 1");
    }
    let query = read_query(accession_query)?;
    let query_values = query
        .replace(',', " ")
        .split_whitespace()
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let requested = parse_accessions(&query_values, version)
        .into_iter()
        .filter(|value| value != "NA")
        .collect::<HashSet<_>>();
    if requested.is_empty() {
        bail!("No accessions were found in --accessions");
    }
    let output_path = output_path.as_ref();
    if let Some(parent) = output_path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    let mut output = BufWriter::with_capacity(IO_BUFFER_BYTES, File::create(output_path)?);
    let mut stats = GrepStats {
        requested: requested.len(),
        ..Default::default()
    };
    let mut records = FastaReader::new(BufReader::with_capacity(
        IO_BUFFER_BYTES,
        File::open(input_path)?,
    ));
    loop {
        cancellation.check()?;
        let batch = record_batch(&mut records, usize::MAX, batch_size)?;
        if batch.is_empty() {
            break;
        }
        grep_batch(
            &mut output,
            &batch,
            &requested,
            version,
            verbose,
            &mut stats,
        )?;
    }
    output.flush()?;
    Ok(stats)
}

#[cfg(test)]
fn extend_accession_set(headers: &mut Vec<Vec<u8>>, accessions: &mut HashSet<String>) {
    let parsed = headers
        .par_iter()
        .map(|header| parse_accession(&String::from_utf8_lossy(header), true))
        .filter(|accession| accession != "NA")
        .collect::<HashSet<_>>();
    accessions.extend(parsed);
    headers.clear();
}

#[cfg(test)]
fn collect_filter_accessions(input_path: &Path, batch_size: usize) -> Result<HashSet<String>> {
    let mut reader = BufReader::with_capacity(IO_BUFFER_BYTES, File::open(input_path)?);
    let mut headers = Vec::with_capacity(batch_size);
    let mut accessions = HashSet::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if !line.starts_with(b">") {
            continue;
        }
        headers.push(std::mem::take(&mut line));
        if headers.len() == batch_size {
            extend_accession_set(&mut headers, &mut accessions);
        }
    }
    extend_accession_set(&mut headers, &mut accessions);
    Ok(accessions)
}

#[cfg(test)]
fn filter_batch(
    output: &mut impl Write,
    batch: &[FastaRecord],
    a2t: &HashMap<String, i64>,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    verbose: bool,
    totals: &mut FilterStats,
) -> Result<()> {
    let accessions = batch
        .par_iter()
        .map(|record| parse_accession(&String::from_utf8_lossy(record.header()), true))
        .collect::<Vec<_>>();
    filter_parsed_batch(
        output,
        batch,
        &accessions,
        a2t,
        filter_taxa,
        mode,
        verbose,
        totals,
    )
}

#[allow(clippy::too_many_arguments)]
fn filter_parsed_batch(
    output: &mut impl Write,
    batch: &[FastaRecord],
    accessions: &[String],
    a2t: &HashMap<String, i64>,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    verbose: bool,
    totals: &mut FilterStats,
) -> Result<()> {
    let decisions = accessions
        .par_iter()
        .map(|accession| {
            if accession == "NA" {
                (mode == FilterMode::Remove, true, false)
            } else if let Some(taxid) = a2t.get(accession) {
                let keep = match mode {
                    FilterMode::Keep => filter_taxa.contains(taxid),
                    FilterMode::Remove => !filter_taxa.contains(taxid),
                };
                (keep, false, false)
            } else {
                (mode == FilterMode::Remove, false, true)
            }
        })
        .collect::<Vec<_>>();

    for ((record, accession), (keep, missing_accession, missing_taxid)) in
        batch.iter().zip(accessions).zip(decisions)
    {
        if missing_accession {
            totals.missing_accession += 1;
            if verbose {
                println!(
                    "No accession found: {}",
                    String::from_utf8_lossy(record.header()).trim()
                );
            }
        }
        if missing_taxid {
            totals.missing_taxid += 1;
            if verbose {
                println!(
                    "No taxid found for {accession}: {}",
                    String::from_utf8_lossy(record.header()).trim()
                );
            }
        }
        if keep {
            record.write_to(output)?;
            totals.kept += 1;
        } else {
            totals.removed += 1;
        }
    }
    Ok(())
}

pub fn filter_fasta(
    input_path: impl AsRef<Path>,
    output_path: Option<&Path>,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    batch_size: usize,
    verbose: bool,
) -> Result<FilterStats> {
    let options = TaxutilsOptions::default();
    filter_fasta_with_options(
        input_path,
        output_path,
        filter_taxa,
        mode,
        batch_size,
        verbose,
        &options.save_folder,
        options.wgs,
    )
}

/// Filter a FASTA using an explicit resource directory and WGS lookup policy.
///
/// This form is intended for language bindings and applications that already
/// resolved their configuration and should not reread process environment.
#[allow(clippy::too_many_arguments)]
pub fn filter_fasta_with_options(
    input_path: impl AsRef<Path>,
    output_path: Option<&Path>,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    batch_size: usize,
    verbose: bool,
    save_folder: impl AsRef<Path>,
    wgs: bool,
) -> Result<FilterStats> {
    filter_fasta_with_options_and_cancel(
        input_path,
        output_path,
        filter_taxa,
        mode,
        batch_size,
        verbose,
        save_folder,
        wgs,
        &CancellationToken::default(),
    )
}

#[allow(clippy::too_many_arguments)]
pub fn filter_fasta_with_options_and_cancel(
    input_path: impl AsRef<Path>,
    output_path: Option<&Path>,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    batch_size: usize,
    verbose: bool,
    save_folder: impl AsRef<Path>,
    wgs: bool,
    cancellation: &CancellationToken,
) -> Result<FilterStats> {
    if batch_size < 1 {
        bail!("--batch-size must be at least 1");
    }
    let input_path = input_path.as_ref();
    let destination = output_path.unwrap_or(input_path);
    let index = AccessionTaxidIndex::open(save_folder, wgs, true, cancellation)?;
    write_filtered_fasta_bounded(
        input_path,
        destination,
        index,
        filter_taxa,
        mode,
        batch_size,
        verbose,
        cancellation,
    )
}

#[allow(clippy::too_many_arguments)]
fn write_filtered_fasta_bounded(
    input_path: &Path,
    destination: &Path,
    mut index: AccessionTaxidIndex,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    batch_size: usize,
    verbose: bool,
    cancellation: &CancellationToken,
) -> Result<FilterStats> {
    let output_dir = destination.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_dir)?;
    let temporary = tempfile::Builder::new()
        .prefix(
            destination
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("tu-filter"),
        )
        .suffix(".tmp")
        .rand_bytes(0)
        .tempfile_in(output_dir)?;
    let mut output = BufWriter::with_capacity(IO_BUFFER_BYTES, temporary);
    let mut records = FastaReader::new(BufReader::with_capacity(
        IO_BUFFER_BYTES,
        File::open(input_path)?,
    ));
    let mut totals = FilterStats::default();
    loop {
        cancellation.check()?;
        let batch = record_batch(&mut records, batch_size, usize::MAX)?;
        if batch.is_empty() {
            break;
        }
        let accessions = batch
            .par_iter()
            .map(|record| parse_accession(&String::from_utf8_lossy(record.header()), true))
            .collect::<Vec<_>>();
        let requested = accessions
            .iter()
            .filter(|accession| accession.as_str() != "NA")
            .cloned()
            .collect::<HashSet<_>>();
        let a2t = index.lookup(requested)?;
        filter_parsed_batch(
            &mut output,
            &batch,
            &accessions,
            &a2t,
            filter_taxa,
            mode,
            verbose,
            &mut totals,
        )?;
    }
    output.flush().context("failed to finish filtered FASTA")?;
    let temporary = output.into_inner().map_err(|error| error.into_error())?;
    temporary
        .persist(destination)
        .map_err(|error| error.error)?;
    Ok(totals)
}

#[cfg(test)]
fn write_filtered_fasta(
    input_path: &Path,
    destination: &Path,
    a2t: &HashMap<String, i64>,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    batch_size: usize,
    verbose: bool,
) -> Result<FilterStats> {
    let output_dir = destination.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_dir)?;
    let temporary = tempfile::Builder::new()
        .prefix(
            destination
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("tu-filter"),
        )
        .suffix(".tmp")
        .rand_bytes(0)
        .tempfile_in(output_dir)?;
    let mut output = BufWriter::with_capacity(IO_BUFFER_BYTES, temporary);
    let mut records = FastaReader::new(BufReader::with_capacity(
        IO_BUFFER_BYTES,
        File::open(input_path)?,
    ));
    let mut totals = FilterStats::default();
    loop {
        let batch = record_batch(&mut records, batch_size, usize::MAX)?;
        if batch.is_empty() {
            break;
        }
        filter_batch(
            &mut output,
            &batch,
            a2t,
            filter_taxa,
            mode,
            verbose,
            &mut totals,
        )?;
    }
    output.flush().context("failed to finish filtered FASTA")?;
    let temporary = output.into_inner().map_err(|error| error.into_error())?;
    temporary
        .persist(destination)
        .map_err(|error| error.error)?;
    Ok(totals)
}

pub fn parse_taxa(value: &str, option_name: &str) -> Result<HashSet<i64>> {
    let text = read_query(value)?;
    let mut taxa = HashSet::new();
    for token in text.replace(',', " ").split_whitespace() {
        taxa.insert(
            token
                .parse()
                .with_context(|| format!("Invalid taxid in {option_name}: {token}"))?,
        );
    }
    if taxa.is_empty() {
        bail!("{option_name} did not contain any taxids");
    }
    Ok(taxa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extract_and_grep_preserve_records_across_small_batches() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.fa");
        let accessions = dir.path().join("accs.txt");
        let output = dir.path().join("out.fa");
        fs::write(
            &input,
            ">NC_045512.2 description\nACGT\n>AB12345.1 other\nTT",
        )
        .unwrap();
        assert_eq!(extract_accessions(&input, &accessions, 1).unwrap(), 2);
        assert_eq!(
            fs::read_to_string(&accessions).unwrap(),
            "NC_045512.2\nAB12345.1\n"
        );
        let stats = grep_fasta(&input, "AB12345.1", &output, true, 1, false).unwrap();
        assert_eq!(stats.matched, 1);
        assert_eq!(fs::read_to_string(output).unwrap(), ">AB12345.1 other\nTT");
    }

    #[test]
    fn cancelled_extract_does_not_install_partial_output() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("cancel-in.fa");
        let output = dir.path().join("cancel-out.txt");
        fs::write(&input, ">NC_045512.2 description\nACGT\n").unwrap();
        let cancellation = CancellationToken::default();
        cancellation.cancel();

        let error = extract_accessions_with_cancel(&input, &output, 1, &cancellation)
            .expect_err("cancelled extraction must fail");

        assert!(error.to_string().contains("operation cancelled"));
        assert!(!output.exists());
    }

    #[test]
    fn clean_preserves_non_header_bytes_and_supports_in_place() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("clean-in.fa");
        let output = dir.path().join("clean-out.fa");
        fs::write(
            &input,
            b"preamble\n>nc_045512.2 description\r\nACGT\r\n>AB12345.1 other\nTT",
        )
        .unwrap();
        let expected = b"preamble\n>NC_045512.2\nACGT\r\n>AB12345.1\nTT";

        clean_fasta_headers(&input, Some(output.as_path()), false).unwrap();
        assert_eq!(fs::read(&output).unwrap(), expected);

        clean_fasta_headers(&input, None, false).unwrap();
        assert_eq!(fs::read(&input).unwrap(), expected);
    }

    #[test]
    fn reader_is_bounded_and_preserves_order() {
        let input = b">A00001.1\nA\n>A00002.1\nBB\n>A00003.1\nCCC\n";
        let mut records = FastaReader::new(BufReader::new(&input[..]));
        let first = record_batch(&mut records, 2, usize::MAX).unwrap();
        assert_eq!(first.len(), 2);
        assert!(String::from_utf8_lossy(first[0].header()).contains("A00001.1"));
        let second = record_batch(&mut records, 2, usize::MAX).unwrap();
        assert_eq!(second.len(), 1);
        assert!(String::from_utf8_lossy(second[0].header()).contains("A00003.1"));
    }

    #[test]
    fn filter_accession_scan_deduplicates_and_ignores_unparseable_headers() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("filter-in.fa");
        fs::write(
            &input,
            ">NC_000001.1 first\nA\n>nc_000001.1 duplicate\nC\n>missing_accession\nG\n>AB12345.2\nT\n",
        )
        .unwrap();

        let accessions = collect_filter_accessions(&input, 2).unwrap();

        assert_eq!(
            accessions,
            HashSet::from(["NC_000001.1".to_owned(), "AB12345.2".to_owned()])
        );
    }

    #[test]
    fn prefetched_filter_map_preserves_keep_and_remove_semantics() {
        fn record(text: &str) -> FastaRecord {
            let data = text.as_bytes().to_vec();
            let header_len = data.iter().position(|byte| *byte == b'\n').unwrap() + 1;
            FastaRecord { data, header_len }
        }

        let batch = vec![
            record(">NC_000001.1 target\nA\n"),
            record(">NC_000002.1 other\nC\n"),
            record(">NC_000003.1 unmapped\nG\n"),
            record(">missing_accession\nT\n"),
        ];
        let a2t = HashMap::from([
            ("NC_000001.1".to_owned(), 13),
            ("NC_000002.1".to_owned(), 15),
        ]);

        let mut keep_output = Vec::new();
        let mut keep_stats = FilterStats::default();
        filter_batch(
            &mut keep_output,
            &batch,
            &a2t,
            &HashSet::from([13]),
            FilterMode::Keep,
            false,
            &mut keep_stats,
        )
        .unwrap();
        assert_eq!(keep_output, batch[0].data);
        assert_eq!(
            keep_stats,
            FilterStats {
                kept: 1,
                removed: 3,
                missing_accession: 1,
                missing_taxid: 1,
            }
        );

        let mut remove_output = Vec::new();
        let mut remove_stats = FilterStats::default();
        filter_batch(
            &mut remove_output,
            &batch,
            &a2t,
            &HashSet::from([15]),
            FilterMode::Remove,
            false,
            &mut remove_stats,
        )
        .unwrap();
        assert_eq!(
            remove_output,
            [
                batch[0].data.as_slice(),
                batch[2].data.as_slice(),
                batch[3].data.as_slice(),
            ]
            .concat()
        );
        assert_eq!(
            remove_stats,
            FilterStats {
                kept: 3,
                removed: 1,
                missing_accession: 1,
                missing_taxid: 1,
            }
        );
    }

    #[test]
    fn filter_supports_atomic_in_place_output() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("same.fa");
        let original = b">NC_000001.1 target\nA\n>NC_000002.1 other\nC\n";
        fs::write(&input, original).unwrap();

        let stats = write_filtered_fasta(
            &input,
            &input,
            &HashMap::from([
                ("NC_000001.1".to_owned(), 13),
                ("NC_000002.1".to_owned(), 15),
            ]),
            &HashSet::from([13]),
            FilterMode::Keep,
            10,
            false,
        )
        .unwrap();

        assert_eq!(stats.kept, 1);
        assert_eq!(stats.removed, 1);
        assert_eq!(fs::read(&input).unwrap(), b">NC_000001.1 target\nA\n");
        assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 1);
    }

    #[test]
    fn parallel_grep_is_deterministic() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("many.fa");
        let one = dir.path().join("one.fa");
        let four = dir.path().join("four.fa");
        let contents = (0..200)
            .map(|index| format!(">NC_{index:06}.1 record {index}\nACGT{index}\n"))
            .collect::<String>();
        fs::write(&input, contents).unwrap();
        let query = (0..200)
            .filter(|index| index % 3 == 0)
            .map(|index| format!("NC_{index:06}.1"))
            .collect::<Vec<_>>()
            .join(",");
        let one_stats = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .unwrap()
            .install(|| grep_fasta(&input, &query, &one, true, 97, false))
            .unwrap();
        let four_stats = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .unwrap()
            .install(|| grep_fasta(&input, &query, &four, true, 97, false))
            .unwrap();
        assert_eq!(one_stats, four_stats);
        assert_eq!(fs::read(one).unwrap(), fs::read(four).unwrap());
    }
}
