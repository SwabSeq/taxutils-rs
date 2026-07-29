//! Rust port of the Python `taxutils` package.
//!
//! The main entry point is [`TaxonomicUtils::new`]. Resource locations and
//! defaults match the Python package: `TAXUTILS_GLOBALS`, or `./taxutils/`.

mod accession;
mod resources;
mod taxonomy;

pub use resources::{TaxutilsBuilder, TaxutilsOptions};
pub use taxonomy::{TaxonId, TaxonNode, TaxonomicUtils, TopologyProfile, TopologyStat};

/// Python-compatible convenience constructor.
pub fn taxutils() -> anyhow::Result<TaxonomicUtils> {
    TaxonomicUtils::new(TaxutilsOptions::default())
}
