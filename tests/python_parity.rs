use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::Command;

use flate2::Compression;
use flate2::write::GzEncoder;
use serde_json::{Value, json};
use taxutils::{TaxonomicUtils, TaxutilsOptions, TopologyStat, parse_accessions};

fn python_root() -> PathBuf {
    std::env::var_os("TAXUTILS_PYTHON_ROOT")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .expect("Rust project has a parent directory")
                .join("taxutils")
        })
}

fn prepare_fixture() -> tempfile::TempDir {
    let fixture = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ncbi");
    let temporary = tempfile::tempdir().expect("create fixture directory");
    for filename in ["names.dmp", "nodes.dmp", "targets.json"] {
        fs::copy(fixture.join(filename), temporary.path().join(filename))
            .expect("copy fixture resource");
    }

    let mut source =
        File::open(fixture.join("nucl_gb.accession2taxid.tsv")).expect("open accession fixture");
    let output =
        File::create(temporary.path().join("nucl_gb.accession2taxid.gz")).expect("create gzip");
    let mut encoder = GzEncoder::new(output, Compression::fast());
    std::io::copy(&mut source, &mut encoder).expect("compress accession fixture");
    encoder.finish().expect("finish accession gzip");
    let mut source = File::open(fixture.join("nucl_wgs.accession2taxid.tsv"))
        .expect("open WGS accession fixture");
    let output = File::create(temporary.path().join("nucl_wgs.accession2taxid.gz"))
        .expect("create WGS gzip");
    let mut encoder = GzEncoder::new(output, Compression::fast());
    std::io::copy(&mut source, &mut encoder).expect("compress WGS fixture");
    encoder.finish().expect("finish WGS accession gzip");
    temporary
}

fn run_python(fixture: &Path) -> Value {
    let root = python_root();
    let source = root.join("src");
    assert!(
        source.join("taxutils/taxutils.py").is_file(),
        "Python taxutils not found at {}. Set TAXUTILS_PYTHON_ROOT.",
        root.display()
    );
    let output = Command::new(std::env::var_os("PYTHON").unwrap_or_else(|| "python3".into()))
        .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/python_reference.py"))
        .arg(fixture)
        .env("TAXUTILS_GLOBALS", fixture)
        .env("PYTHONPATH", source)
        .output()
        .expect("execute Python reference");
    assert!(
        output.status.success(),
        "Python reference failed:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
    serde_json::from_slice(&output.stdout).expect("Python emitted valid JSON")
}

fn sorted_a2t(tu: &TaxonomicUtils) -> Vec<(String, i64)> {
    let mut values = tu
        .a2t
        .iter()
        .map(|(accession, taxid)| (accession.clone(), *taxid))
        .collect::<Vec<_>>();
    values.sort();
    values
}

fn rust_results(fixture: &Path) -> Value {
    let mut tu = TaxonomicUtils::new(TaxutilsOptions {
        save_folder: fixture.to_owned(),
        low_memory: true,
        ..Default::default()
    })
    .expect("load Rust fixture");

    let accession_inputs = [
        ">NC_000001.1 bacterial alpha",
        ">kraken:taxid|22|ABCD01000001.1 viral example",
        "prefix NC_000003.2 suffix",
        "nc_000002.1",
        "sequence_without_accession",
        "FOONC_045512.2",
    ];
    tu.load_a2t(
        &[">NC_000001.1 description", "ABCD01000001.1", "missing"],
        None,
        false,
        None,
    )
    .expect("replace accession map");
    let a2t_replace = sorted_a2t(&tu);
    tu.load_a2t(&["NC_000002.1", "NC_000001.1"], None, true, None)
        .expect("extend accession map");
    let a2t_extend = sorted_a2t(&tu);
    tu.load_a2t(&["WXYZ01000001.1"], None, false, Some(true))
        .expect("load WGS accession map");
    let a2t_wgs = sorted_a2t(&tu);
    tu.load_a2t(
        &["NC_000001.1", "WXYZ01000001.1"],
        Some(false),
        false,
        Some(true),
    )
    .expect("load SQLite WGS accession map");
    let a2t_sqlite_wgs = sorted_a2t(&tu);

    let topology_many = tu.topologies(&[13, 15], Some("F")).expect("batch topology");
    let mut t2a = tu
        .get_t2a(&[13, 22], None, None)
        .expect("taxid to accession")
        .into_iter()
        .collect::<Vec<_>>();
    t2a.sort();
    let mut t2a_sqlite_wgs = tu
        .get_t2a(&[15, 22], Some(false), Some(true))
        .expect("SQLite WGS taxid to accession")
        .into_iter()
        .collect::<Vec<_>>();
    t2a_sqlite_wgs.sort();
    let topology_stats = json!({
        "n_taxa": tu.topology_stat(12, None, TopologyStat::NTaxa).unwrap(),
        "n_leaves": tu.topology_stat(12, None, TopologyStat::NLeaves).unwrap(),
        "max_depth": tu.topology_stat(12, None, TopologyStat::MaxDepth).unwrap(),
        "mean_depth": tu.topology_stat(12, None, TopologyStat::MeanDepth).unwrap(),
        "topology_scale": tu.topology_stat(12, None, TopologyStat::TopologyScale).unwrap(),
        "max_children": tu.topology_stat(12, None, TopologyStat::MaxChildren).unwrap(),
        "branching_taxa_fraction": tu.topology_stat(
            12, None, TopologyStat::BranchingTaxaFraction
        ).unwrap(),
        "top_child_fraction": tu.topology_stat(
            12, None, TopologyStat::TopChildFraction
        ).unwrap(),
    });

    let node_taxa = [1, 2, 3, 10, 11, 12, 13, 14, 15, 10239, 20, 21, 22];
    let names = [0, 1, 2, 13, 22, 694009, 2697049].map(|taxon| json!([taxon, tu.names[&taxon]]));
    let parents = node_taxa.map(|taxon| json!([taxon, tu.parent[&taxon]]));
    let leaf_many = [12, 14, 15, 999999].map(|taxon| tu.is_leaf(taxon));
    let child_many =
        [(14, 13), (15, 12), (22, 21)].map(|(child, parent)| tu.is_child(child, parent));
    let descendent_many = [(14, 10), (13, 13), (10, 10), (22, 20)]
        .map(|(child, parent)| tu.is_descendent(child, parent));

    json!({
        "accessions_versioned": parse_accessions(&accession_inputs, true),
        "accessions_unversioned": parse_accessions(&accession_inputs, false),
        "rank_order": tu.get_rank_order(),
        "names": names,
        "nodes": tu.nodes,
        "parents": parents,
        "target_taxa": tu.target_taxa,
        "branch": tu.get_branch(14),
        "subtree": tu.get_subtree(12),
        "ancestors": tu.get_ancestors(&[14, 15, 22], "F").unwrap(),
        "leaf_scalar": tu.is_leaf(14),
        "leaf_many": leaf_many,
        "child_scalar": tu.is_child(14, 13),
        "child_many": child_many,
        "descendent_scalar": tu.is_descendent(14, 10),
        "descendent_many": descendent_many,
        "lca": tu.get_lca(14, 15),
        "distance": tu.get_distance(14, 15),
        "sort_taxa": tu.sort_taxa([22, 15, 10, 14, 13, 15]),
        "tree_with_ancestors": tu.format_tree([14, 15, 22], true, 1, "\t"),
        "tree_without_ancestors": tu.format_tree([14, 15, 22], false, 1, "\t"),
        "higher_than_family": tu.higher_than_rank(
            &[1, 2, 10, 11, 12, 13, 999999],
            "F"
        ).unwrap(),
        "topology": tu.topology(12, None).unwrap(),
        "topology_family_anchor": tu.topology(14, Some("F")).unwrap(),
        "topology_many": topology_many,
        "topology_scale": topology_many
            .iter()
            .map(|profile| profile.topology_scale)
            .collect::<Vec<_>>(),
        "topology_stats": topology_stats,
        "a2t_replace": a2t_replace,
        "a2t_extend": a2t_extend,
        "a2t_wgs": a2t_wgs,
        "a2t_sqlite_wgs": a2t_sqlite_wgs,
        "t2a": t2a,
        "t2a_sqlite_wgs": t2a_sqlite_wgs,
    })
}

#[test]
fn core_api_matches_python_taxutils() {
    let fixture = prepare_fixture();
    let python = run_python(fixture.path());
    let rust = rust_results(fixture.path());
    assert_eq!(
        rust,
        python,
        "Rust/Python parity mismatch.\nRust:\n{}\nPython:\n{}",
        serde_json::to_string_pretty(&rust).unwrap(),
        serde_json::to_string_pretty(&python).unwrap()
    );
}
