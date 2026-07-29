# taxutils-rs

Rust port of Python [`taxutils`](../taxutils), preserving its NCBI taxonomy,
accession lookup, corrected-rank, target-taxa, topology, and FASTA functionality.

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
otherwise uses `./taxutils/`, matching the Python package. It manages:

- `names.dmp` and `nodes.dmp` from the NCBI taxdump
- `targets.json`
- `nucl_gb.accession2taxid.gz`
- optional `nucl_wgs.accession2taxid.gz`
- `nucl.accession2taxid.db` for `low_memory = false`

The public API includes accession parsing and bidirectional accession/taxid
lookups; branches, subtrees, ancestors, leaves, child/descendant checks, LCAs,
distances, taxonomic sorting and tree formatting; corrected ranks and rank
threshold checks; target taxa; and all topology metrics from Python 1.0.3.

Rust returns typed values (`Vec`, `HashMap`, `HashSet`, `TaxonNode`, and
`TopologyProfile`) where Python returns pandas or NumPy containers.

## Command line

```console
tu extract input.fasta -o accessions.txt
tu clean -i input.fasta [-o clean.fasta]
tu grep -i input.fasta -a NC_045512.2,NC_001422.1 -o hits.fasta
tu filter -i input.fasta -o filtered.fasta --keep-taxids 2697049
```

Taxid and accession arguments may also name text files. `grep --no-version`
matches accessions without versions. `filter` uses the indexed SQLite mode just
like the Python CLI.

## Verification

```console
cargo test --all-targets
cargo clippy --all-targets --all-features -- -D warnings
```

The tests cover the Python accession examples and boundaries, corrected rank
assignment, tree relationships, LCA/distance, topology formulae, and FASTA
record behavior.
