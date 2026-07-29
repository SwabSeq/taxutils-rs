use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};
use rayon::prelude::*;

use crate::{TaxonomicUtils, TaxutilsOptions, parse_accession, parse_accessions};

const CLEAN_BATCH_BYTES: usize = 1_000_000;

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
            let mut line = Vec::new();
            match self.reader.read_until(b'\n', &mut line) {
                Ok(0) => {
                    self.finished = true;
                    break;
                }
                Ok(_) if line.starts_with(b">") => {
                    self.pending_header = Some(line);
                    break;
                }
                Ok(_) => data.extend_from_slice(&line),
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
    headers: &mut Vec<(usize, String)>,
) -> Result<usize> {
    let accessions = parse_accessions(
        &headers.iter().map(|(_, header)| header).collect::<Vec<_>>(),
        true,
    );
    for ((line_number, header), accession) in headers.iter().zip(&accessions) {
        if accession == "NA" {
            bail!(
                "No accession found in FASTA header on line {line_number}: {}",
                header.trim()
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
    if batch_size < 1 {
        bail!("--batch-size must be at least 1");
    }
    let output_path = output_path.as_ref();
    let output_dir = output_path.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_dir)?;
    let mut temporary = tempfile::NamedTempFile::new_in(output_dir)?;
    let mut headers = Vec::with_capacity(batch_size);
    let mut count = 0;
    for (line_number, line) in BufReader::new(File::open(fasta_path)?).lines().enumerate() {
        let line = line?;
        if !line.starts_with('>') {
            continue;
        }
        headers.push((line_number + 1, line));
        if headers.len() == batch_size {
            count += write_accession_batch(&mut temporary, &mut headers)?;
        }
    }
    count += write_accession_batch(&mut temporary, &mut headers)?;
    temporary.flush()?;
    temporary
        .persist(output_path)
        .map_err(|error| error.error)?;
    Ok(count)
}

fn write_clean_batch(
    output: &mut impl Write,
    lines: &mut Vec<(usize, Vec<u8>)>,
    verbose: bool,
) -> Result<()> {
    let accessions = lines
        .par_iter()
        .map(|(_, line)| {
            line.starts_with(b">")
                .then(|| parse_accession(&String::from_utf8_lossy(line), true))
        })
        .collect::<Vec<_>>();
    for ((line_number, line), accession) in lines.iter().zip(accessions) {
        if let Some(accession) = accession {
            if accession == "NA" {
                if verbose {
                    println!(
                        "NA accession line {line_number}: {}",
                        String::from_utf8_lossy(line).trim()
                    );
                }
                bail!(
                    "No accession found in FASTA header on line {line_number}: {}",
                    String::from_utf8_lossy(line).trim()
                );
            }
            writeln!(output, ">{accession}")?;
        } else {
            output.write_all(line)?;
        }
    }
    lines.clear();
    Ok(())
}

pub fn clean_fasta_headers(
    input_path: impl AsRef<Path>,
    output_path: Option<&Path>,
    verbose: bool,
) -> Result<()> {
    let input_path = input_path.as_ref();
    let destination = output_path.unwrap_or(input_path);
    let output_dir = destination.parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(output_dir)?;
    let mut temporary = tempfile::NamedTempFile::new_in(output_dir)?;
    let mut reader = BufReader::new(File::open(input_path)?);
    let mut lines = Vec::new();
    let mut buffered_bytes = 0;
    let mut line_number = 0;
    loop {
        let mut line = Vec::new();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        line_number += 1;
        buffered_bytes += line.len();
        lines.push((line_number, line));
        if buffered_bytes >= CLEAN_BATCH_BYTES {
            write_clean_batch(&mut temporary, &mut lines, verbose)?;
            buffered_bytes = 0;
        }
    }
    write_clean_batch(&mut temporary, &mut lines, verbose)?;
    temporary.flush()?;
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
    let mut output = BufWriter::new(File::create(output_path)?);
    let mut stats = GrepStats {
        requested: requested.len(),
        ..Default::default()
    };
    let mut records = FastaReader::new(BufReader::new(File::open(input_path)?));
    loop {
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

fn filter_batch(
    output: &mut impl Write,
    batch: &[FastaRecord],
    tu: &mut TaxonomicUtils,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    verbose: bool,
    totals: &mut FilterStats,
) -> Result<()> {
    let accessions = batch
        .par_iter()
        .map(|record| parse_accession(&String::from_utf8_lossy(record.header()), true))
        .collect::<Vec<_>>();
    let lookup = accessions
        .iter()
        .filter(|accession| *accession != "NA")
        .collect::<HashSet<_>>()
        .into_iter()
        .cloned()
        .collect::<Vec<_>>();
    if !lookup.is_empty() {
        tu.load_a2t(&lookup, None, false, None)?;
    }
    let decisions = accessions
        .par_iter()
        .map(|accession| {
            if accession == "NA" {
                (mode == FilterMode::Remove, true, false)
            } else if let Some(taxid) = tu.a2t.get(accession) {
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
    output_path: impl AsRef<Path>,
    filter_taxa: &HashSet<i64>,
    mode: FilterMode,
    batch_size: usize,
    verbose: bool,
) -> Result<FilterStats> {
    if batch_size < 1 {
        bail!("--batch-size must be at least 1");
    }
    let output_path = output_path.as_ref();
    if let Some(parent) = output_path
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    let mut tu = TaxonomicUtils::new(TaxutilsOptions {
        low_memory: false,
        ..Default::default()
    })?;
    let mut output = BufWriter::new(File::create(output_path)?);
    let mut records = FastaReader::new(BufReader::new(File::open(input_path)?));
    let mut totals = FilterStats::default();
    loop {
        let batch = record_batch(&mut records, batch_size, usize::MAX)?;
        if batch.is_empty() {
            break;
        }
        filter_batch(
            &mut output,
            &batch,
            &mut tu,
            filter_taxa,
            mode,
            verbose,
            &mut totals,
        )?;
    }
    output.flush().context("failed to finish filtered FASTA")?;
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
