# taxutils-rs

Rust port of the Python
[`taxutils`](https://github.com/SwabSeq/taxutils) package, optimized for parallelism and incorporation into Rust repositories with its [crate](https://crates.io/crates/taxutils). The port was performed entirely by codex with GPT-5.6-sol. Available for download from [bioconda](https://anaconda.org/channels/bioconda/packages/taxutils-rs/overview).

The Rust package provides both:

- the `tu` command-line program
- the `taxutils` Rust library crate, [crates.io]

## Installation

Install via conda or mamba from the bioconda forge:
```console
conda install bioconda::taxutils-rs
```

Or install the command-line program from crates.io:

```console
cargo install taxutils
tu --help
```

The Cargo package and library are named `taxutils`; the installed executable is
named `tu`. To use the library in another Rust project:

```console
cargo add taxutils
```

## Data directory

Like the Python package, taxutils-rs stores downloaded NCBI taxonomy and
accession resources in `./taxutils/` by default. Set `TAXUTILS_GLOBALS` to use a persistent location instead:

```console
export TAXUTILS_GLOBALS=/path/to/taxutils/cache
tu filter -i input.fasta -o filtered.fasta --keep-taxids 2697049
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
use taxutils::{TaxonomicUtils, TaxutilsBuilder, TopologyStat};

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
let scale = tu.topology_stat(
    2697049, Some("F"), TopologyStat::TopologyScale)?;
# Ok::<(), anyhow::Error>(())
```

`TaxutilsOptions::default()` reads `TAXUTILS_GLOBALS` at construction time and
otherwise uses `./taxutils/`, matching the Python package.

The library API is centered on `TaxonomicUtils` and provides the Python object
capabilities:

- `parse_accession`, `load_a2t`, and `get_t2a`
- `get_branch`, `get_subtree`, `get_ancestor`, `is_leaf`, `is_child`, and
  `is_descendent`
- `get_lca`, `get_distance`, `sort_taxa`, and `format_tree`
- `topology` and `topology_stat`
- `get_rank_order` and `higher_than_rank`

Python uses one method name for both scalar and pandas/NumPy inputs. Rust keeps
the same scalar names and provides explicit, parallel batch forms:

| Python capability | Rust scalar | Rust batch |
| --- | --- | --- |
| `parse_accession` | `parse_accession` | `parse_accessions` |
| `get_ancestor` | `get_ancestor` | `get_ancestors` |
| `is_leaf` | `is_leaf` | `are_leaves` |
| `is_child` | `is_child` | `are_children` |
| `is_descendent` | `is_descendent` | `are_descendents` |
| `topology` | `topology` | `topologies` |
| `topology(..., stat=...)` | `topology_stat` | `topology_stats` |

The object exposes `names`, `nodes`, `parent`, `target_taxa`, and `a2t`, using
typed Rust collections and records (`Vec`, `HashMap`, `HashSet`, `TaxonNode`,
and `TopologyProfile`) where Python returns pandas or NumPy containers.

FASTA extraction, cleaning, grepping, and filtering are intentionally not
library functions. They are available only through the `tu` commands below.

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
to cap worker threads. FASTA records are written in bounded, ordered batches,
so parallel execution does not reorder records. `filter` first collects one
entry per unique accession for a single bulk SQLite lookup, then scans the FASTA
again to write records.

Taxid and accession arguments may also name text files. `grep --no-version`
matches accessions without versions. `filter` uses the indexed SQLite mode just
like the Python CLI and is the command that uses the resources under
`TAXUTILS_GLOBALS`. `extract`, `clean`, and `grep` operate directly on FASTA
data and do not require the NCBI resource cache.

## Performance model

- FASTA records are streamed in bounded batches instead of retaining sequence
  data in memory. `filter` retains one map entry per unique parsed accession.
- Header parsing and filter decisions run in parallel; large buffered reads and
  writes preserve input order without issuing I/O for every FASTA line.
- `filter` performs one bulk accession-to-taxid query before its output pass,
  avoiding a SQLite connection, temporary table, and join for every batch.
- `clean` is normally limited by sequential input/output throughput. Additional
  threads accelerate header parsing but do not divide the input file into shards.
- `--batch-size` controls records per batch for `extract` and `filter`, and
  approximate bytes per batch for `grep`.
- GB and WGS gzip mapping files are scanned concurrently when WGS mode is enabled.
- SQLite lookups use indexed temporary-table joins instead of one query per key.
- NCBI dump parsing avoids allocating a temporary field vector for every row.
- Batch accession parsing, ancestor lookup, rank checks, target expansion, node
  materialization, and `topologies` use the shared Rayon worker pool.
