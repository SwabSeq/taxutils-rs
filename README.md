# taxutils-rs

Rust port of the Python
[`taxutils`](https://pypi.org/project/taxutils/) package, preserving its NCBI
taxonomy, accession lookup, corrected-rank, target-taxa, topology, and FASTA
functionality. The original Python package is available from
[PyPI](https://pypi.org/project/taxutils/) and
[Bioconda](https://anaconda.org/channels/bioconda/packages/taxutils/overview).

The Rust package provides both:

- the `tu` command-line program
- the `taxutils` Rust library crate

## Installation

Install the command-line program from crates.io:

```console
cargo install taxutils
tu --help
```

The Cargo package and library are named `taxutils`; the installed executable is
named `tu`. To use the library in another Rust project:

```console
cargo add taxutils
```

`taxutils-rs` is also intended to be distributed through Bioconda under the
package name `taxutils-rs`. Until that recipe is published, install from source:

```console
git clone https://github.com/SwabSeq/taxutils-rs.git
cargo install --path taxutils-rs
```

## Data directory

Like the Python package, taxutils-rs stores downloaded NCBI taxonomy and
accession resources in `./taxutils/` by default. Set `TAXUTILS_GLOBALS` to use
a persistent location instead:

```console
export TAXUTILS_GLOBALS=/path/to/taxutils/cache
tu filter -i input.fasta -o filtered.fasta --keep-taxids 2697049
```

The directory is created if needed and should be writable. Keeping it outside
the current working directory prevents repeated downloads when commands are run
from different projects. You can also set it for a single command:

```console
TAXUTILS_GLOBALS=/path/to/taxutils/cache \
  tu filter -i input.fasta -o filtered.fasta --remove-taxids 9606
```

The Rust library reads `TAXUTILS_GLOBALS` when `TaxutilsOptions::default()` or
`TaxutilsBuilder::new()` is constructed. Set the variable before constructing
the builder. Applications that manage configuration directly can set the same
path without an environment variable:

```rust
let tu = taxutils::TaxutilsBuilder::new()
    .save_folder("/path/to/taxutils/cache")
    .build()?;
# Ok::<(), anyhow::Error>(())
```

Managed resources include:

- `names.dmp` and `nodes.dmp` from the NCBI taxdump
- `targets.json`
- `nucl_gb.accession2taxid.gz`
- optional `nucl_wgs.accession2taxid.gz`
- `nucl.accession2taxid.db` when using indexed SQLite lookup

## Library

```rust
use taxutils::{TaxonomicUtils, TaxutilsBuilder};

let mut tu = TaxutilsBuilder::new()
    .low_memory(true)
    .wgs(false)
    .build()?;

let accession = tu.parse_accession(">NC_045512.2 SARS-CoV-2", true);
tu.load_a2t(&[accession], None, false, None)?;

let branch = tu.get_branch(2697049);
let subtree = tu.get_subtree(694009);
let lca = tu.get_lca(2697049, 694009);
let distance = tu.get_distance(2697049, 694009);
let profile = tu.topology(2697049, Some("F"))?;
# Ok::<(), anyhow::Error>(())
```

`TaxutilsOptions::default()` reads `TAXUTILS_GLOBALS` at construction time and
otherwise uses `./taxutils/`, matching the Python package.

The public API includes accession parsing and bidirectional accession/taxid
lookups; branches, subtrees, ancestors, leaves, child/descendant checks, LCAs,
distances, taxonomic sorting and tree formatting; corrected ranks and rank
threshold checks; target taxa; and all topology metrics from Python 1.0.3.

Rust returns typed values (`Vec`, `HashMap`, `HashSet`, `TaxonNode`, and
`TopologyProfile`) where Python returns pandas or NumPy containers.

## Command line

```console
tu extract input.fasta -o accessions.txt
tu clean -i input.fasta -o clean.fasta
tu grep -i input.fasta -a NC_045512.2,NC_001422.1 -o hits.fasta
tu filter -i input.fasta -o filtered.fasta --keep-taxids 2697049
tu filter -i input.fasta -o filtered.fasta --remove-taxids 9606
```

The CLI uses all available logical CPUs for accession parsing and filtering.
Use `tu --threads N <command> ...` (or place `--threads N` after the command)
to cap worker threads. FASTA processing is streamed in bounded, ordered batches,
so parallel execution does not reorder records or load the entire input file.

Taxid and accession arguments may also name text files. `grep --no-version`
matches accessions without versions. `filter` uses the indexed SQLite mode just
like the Python CLI and is the command that uses the resources under
`TAXUTILS_GLOBALS`. `extract`, `clean`, and `grep` operate directly on FASTA
data and do not require the NCBI resource cache.

## Performance model

- FASTA commands stream bounded batches instead of retaining the whole input.
- Header parsing and filter decisions run in parallel; records are written serially
  in input order for deterministic FASTA output.
- `--batch-size` controls records per batch for `extract` and `filter`, and
  approximate bytes per batch for `grep`.
- GB and WGS gzip mapping files are scanned concurrently when WGS mode is enabled.
- SQLite lookups use indexed temporary-table joins instead of one query per key.
- NCBI dump parsing avoids allocating a temporary field vector for every row.
- Batch accession parsing, ancestor lookup, rank checks, target expansion, node
  materialization, and `topologies` use the shared Rayon worker pool.

## Verification

```console
cargo test --all-targets
cargo clippy --all-targets --all-features -- -D warnings
```

### Python parity suite

The integration suite builds a miniature NCBI dataset and compares the Rust
results directly with the neighboring Python `taxutils` package:

```console
tests/run_parity.sh
```

By default it expects Python taxutils at `../taxutils`. Override the package or
interpreter when needed:

```console
TAXUTILS_PYTHON_ROOT=/path/to/taxutils PYTHON=/path/to/python \
    tests/run_parity.sh
```

The API contract covers accession parsing; names, nodes, corrected ranks,
parents, and targets; branches, subtrees, ancestors, leaves, child/descendant
checks, LCA, distance, sorting, tree formatting, rank thresholds, all topology
statistics, low-memory and SQLite accession maps, WGS upgrades, and reverse
taxid lookup. CLI parity covers `extract`, `clean`, versioned and unversioned
`grep`, and both `filter` modes, comparing FASTA files and summary output.
