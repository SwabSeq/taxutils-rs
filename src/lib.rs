//! Rust port of the Python `taxutils` package.
//!
//! The main entry point is [`TaxonomicUtils::new`]. Resource locations and
//! defaults match the Python package: `TAXUTILS_GLOBALS`, or `./taxutils/`.

mod accession;
mod fasta;
mod resources;
mod taxonomy;

pub use accession::{parse_accession, parse_accessions};
pub use fasta::{
    FilterMode, FilterStats, GrepStats, clean_fasta_headers, extract_accessions, filter_fasta,
    grep_fasta, parse_taxa,
};
pub use resources::{TaxutilsBuilder, TaxutilsOptions};
pub use taxonomy::{RANK_CODES, TaxonId, TaxonNode, TaxonomicUtils, TopologyProfile, TopologyStat};

/// Python-compatible convenience constructor.
pub fn taxutils() -> anyhow::Result<TaxonomicUtils> {
    TaxonomicUtils::new(TaxutilsOptions::default())
}
