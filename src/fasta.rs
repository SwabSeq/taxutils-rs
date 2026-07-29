use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::{TaxonomicUtils, TaxutilsOptions, parse_accession};

type FastaRecord = (Vec<u8>, Vec<Vec<u8>>);

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

fn records(path: &Path) -> Result<Vec<FastaRecord>> {
    let mut result = Vec::new();
    let mut header: Option<Vec<u8>> = None;
    let mut sequence = Vec::new();
    let mut reader = BufReader::new(File::open(path)?);
    loop {
        let mut line = Vec::new();
        if reader.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.starts_with(b">") {
            if let Some(previous) = header.replace(line) {
                result.push((previous, std::mem::take(&mut sequence)));
            }
        } else if header.is_some() {
            sequence.push(line);
        }
    }
    if let Some(header) = header {
        result.push((header, sequence));
    }
    Ok(result)
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
    let mut count = 0;
    for (line_number, line) in BufReader::new(File::open(fasta_path)?).lines().enumerate() {
        let line = line?;
        if !line.starts_with('>') {
            continue;
        }
        let accession = parse_accession(&line, true);
        if accession == "NA" {
            bail!(
                "No accession found in FASTA header on line {}: {}",
                line_number + 1,
                line
            );
        }
        writeln!(temporary, "{accession}")?;
        count += 1;
    }
    temporary.flush()?;
    temporary
        .persist(output_path)
        .map_err(|error| error.error)?;
    Ok(count)
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
    for (line_number, line) in BufReader::new(File::open(input_path)?).lines().enumerate() {
        let line = line?;
        if line.starts_with('>') {
            let accession = parse_accession(&line, true);
            if accession == "NA" {
                if verbose {
                    println!("NA accession line {}: {}", line_number + 1, line);
                }
                bail!(
                    "No accession found in FASTA header on line {}: {}",
                    line_number + 1,
                    line
                );
            }
            writeln!(temporary, ">{accession}")?;
        } else {
            writeln!(temporary, "{line}")?;
        }
    }
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

pub fn grep_fasta(
    input_path: impl AsRef<Path>,
    accession_query: &str,
    output_path: impl AsRef<Path>,
    version: bool,
    _batch_size: usize,
    verbose: bool,
) -> Result<GrepStats> {
    let requested = read_query(accession_query)?
        .replace(',', " ")
        .split_whitespace()
        .map(|value| parse_accession(value, version))
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
    for (header, sequence) in records(input_path.as_ref())? {
        stats.scanned += 1;
        let text = String::from_utf8_lossy(&header);
        let accession = parse_accession(&text, version);
        if accession == "NA" {
            stats.missing_accession += 1;
            if verbose {
                println!("NA accession: {}", text.trim());
            }
        } else if requested.contains(&accession) {
            output.write_all(&header)?;
            for line in sequence {
                output.write_all(&line)?;
            }
            stats.matched += 1;
        }
    }
    output.flush()?;
    Ok(stats)
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
    let all_records = records(input_path.as_ref())?;
    let mut output = BufWriter::new(File::create(output_path)?);
    let mut totals = FilterStats::default();
    for chunk in all_records.chunks(batch_size) {
        let accessions = chunk
            .iter()
            .map(|(header, _)| parse_accession(&String::from_utf8_lossy(header), true))
            .collect::<Vec<_>>();
        let lookup = accessions
            .iter()
            .filter(|accession| *accession != "NA")
            .cloned()
            .collect::<Vec<_>>();
        if !lookup.is_empty() {
            tu.load_a2t(&lookup, None, false, None)?;
        }
        for ((header, sequence), accession) in chunk.iter().zip(accessions) {
            let keep = if accession == "NA" {
                totals.missing_accession += 1;
                if verbose {
                    println!(
                        "No accession found: {}",
                        String::from_utf8_lossy(header).trim()
                    );
                }
                mode == FilterMode::Remove
            } else if let Some(taxid) = tu.a2t.get(&accession) {
                match mode {
                    FilterMode::Keep => filter_taxa.contains(taxid),
                    FilterMode::Remove => !filter_taxa.contains(taxid),
                }
            } else {
                totals.missing_taxid += 1;
                if verbose {
                    println!(
                        "No taxid found for {accession}: {}",
                        String::from_utf8_lossy(header).trim()
                    );
                }
                mode == FilterMode::Remove
            };
            if keep {
                output.write_all(header)?;
                for line in sequence {
                    output.write_all(line)?;
                }
                totals.kept += 1;
            } else {
                totals.removed += 1;
            }
        }
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
    fn extract_and_grep_preserve_records() {
        let dir = tempfile::tempdir().unwrap();
        let input = dir.path().join("in.fa");
        let accessions = dir.path().join("accs.txt");
        let output = dir.path().join("out.fa");
        fs::write(
            &input,
            ">NC_045512.2 description\nACGT\n>AB12345.1 other\nTT\n",
        )
        .unwrap();
        assert_eq!(extract_accessions(&input, &accessions, 1).unwrap(), 2);
        assert_eq!(
            fs::read_to_string(&accessions).unwrap(),
            "NC_045512.2\nAB12345.1\n"
        );
        let stats = grep_fasta(&input, "AB12345.1", &output, true, 10, false).unwrap();
        assert_eq!(stats.matched, 1);
        assert_eq!(
            fs::read_to_string(output).unwrap(),
            ">AB12345.1 other\nTT\n"
        );
    }
}
