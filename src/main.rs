use std::path::PathBuf;

use anyhow::Result;
use clap::{Parser, Subcommand};
use taxutils::{
    FilterMode, clean_fasta_headers, extract_accessions, filter_fasta, grep_fasta, parse_taxa,
};

#[derive(Parser)]
#[command(
    name = "tu",
    version,
    about = "Utilities for working with taxonomy data and FASTA files."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Extract one accession per FASTA header.
    Extract {
        fasta: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long, default_value_t = 10_000)]
        batch_size: usize,
    },
    /// Replace FASTA headers with accession-only headers.
    Clean {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(long)]
        verbose: bool,
    },
    /// Extract FASTA records matching requested accessions.
    Grep {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        accessions: String,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(long)]
        no_version: bool,
        #[arg(long, default_value_t = 1_000_000)]
        batch_size: usize,
        #[arg(long)]
        verbose: bool,
    },
    /// Filter FASTA records using accession-to-taxid lookup.
    Filter {
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long)]
        output: PathBuf,
        #[arg(
            long,
            conflicts_with = "remove_taxids",
            required_unless_present = "remove_taxids"
        )]
        keep_taxids: Option<String>,
        #[arg(
            long,
            conflicts_with = "keep_taxids",
            required_unless_present = "keep_taxids"
        )]
        remove_taxids: Option<String>,
        #[arg(long, default_value_t = 5000)]
        batch_size: usize,
        #[arg(long)]
        verbose: bool,
    },
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Extract {
            fasta,
            output,
            batch_size,
        } => {
            let count = extract_accessions(fasta, &output, batch_size)?;
            println!("Wrote {count} accessions to {}", output.display());
        }
        Command::Clean {
            input,
            output,
            verbose,
        } => clean_fasta_headers(input, output.as_deref(), verbose)?,
        Command::Grep {
            input,
            accessions,
            output,
            no_version,
            batch_size,
            verbose,
        } => {
            let totals = grep_fasta(input, &accessions, output, !no_version, batch_size, verbose)?;
            println!(
                "Finished grepping FASTA: requested={} scanned={} matched={} missing_accession={}",
                totals.requested, totals.scanned, totals.matched, totals.missing_accession
            );
        }
        Command::Filter {
            input,
            output,
            keep_taxids,
            remove_taxids,
            batch_size,
            verbose,
        } => {
            let (value, option, mode) = if let Some(value) = keep_taxids {
                (value, "--keep-taxids", FilterMode::Keep)
            } else {
                (
                    remove_taxids.expect("clap requires one filtering option"),
                    "--remove-taxids",
                    FilterMode::Remove,
                )
            };
            let taxa = parse_taxa(&value, option)?;
            let totals = filter_fasta(input, output, &taxa, mode, batch_size, verbose)?;
            println!(
                "Finished filtering FASTA: kept={} removed={} missing_accession={} missing_taxid={}",
                totals.kept, totals.removed, totals.missing_accession, totals.missing_taxid
            );
        }
    }
    Ok(())
}
