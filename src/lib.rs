//! Rust port of the Python `taxutils` package.
//!
//! The main entry point is [`TaxonomicUtils::new`]. Resource locations and
//! defaults match the Python package: `TAXUTILS_GLOBALS`, or `./taxutils/`.

mod accession;
pub mod fasta;
mod resources;
mod taxonomy;
pub mod threads;

pub use fasta::{
    CancellationToken, FilterMode, FilterStats, GrepStats, clean_fasta_headers,
    clean_fasta_headers_with_cancel, extract_accessions, extract_accessions_with_cancel,
    filter_fasta, filter_fasta_with_options, filter_fasta_with_options_and_cancel, grep_fasta,
    grep_fasta_with_cancel, parse_taxa,
};
pub use resources::{
    AccessionDatabaseOptions, TaxutilsBuilder, TaxutilsOptions, ensure_accession_database,
    ensure_accession_database_with_cancel, lookup_accession_taxids,
    lookup_accession_taxids_with_cancel, lookup_taxid_accessions,
    lookup_taxid_accessions_with_cancel,
};
pub use taxonomy::{TaxonId, TaxonNode, TaxonomicUtils, TopologyProfile, TopologyStat};

/// Version of the Rust implementation backing this library.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Python-compatible convenience constructor.
pub fn taxutils() -> anyhow::Result<TaxonomicUtils> {
    TaxonomicUtils::new(TaxutilsOptions::default())
}
