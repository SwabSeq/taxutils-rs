#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
python_root="${TAXUTILS_PYTHON_ROOT:-$(dirname "$repo_root")/taxutils}"
python_bin="${PYTHON:-python3}"
tu_bin="${TU_BIN:-$repo_root/target/debug/tu}"
fixtures="$repo_root/tests/fixtures"

if [[ ! -f "$python_root/src/taxutils/taxutils.py" ]]; then
    echo "Python taxutils not found at $python_root; set TAXUTILS_PYTHON_ROOT" >&2
    exit 1
fi
if [[ ! -x "$tu_bin" ]]; then
    cargo build --manifest-path "$repo_root/Cargo.toml" --bin tu
fi

work="$(mktemp -d "${TMPDIR:-/tmp}/taxutils-parity.XXXXXX")"
trap 'rm -rf "$work"' EXIT

prepare_globals() {
    local destination="$1"
    mkdir -p "$destination"
    cp "$fixtures/ncbi/names.dmp" "$destination/names.dmp"
    cp "$fixtures/ncbi/nodes.dmp" "$destination/nodes.dmp"
    cp "$fixtures/ncbi/targets.json" "$destination/targets.json"
    gzip -c "$fixtures/ncbi/nucl_gb.accession2taxid.tsv" \
        > "$destination/nucl_gb.accession2taxid.gz"
    gzip -c "$fixtures/ncbi/nucl_wgs.accession2taxid.tsv" \
        > "$destination/nucl_wgs.accession2taxid.gz"
}

prepare_globals "$work/python-globals"
prepare_globals "$work/rust-globals"
mkdir -p "$work/python" "$work/rust"

run_python() {
    TAXUTILS_GLOBALS="$work/python-globals" \
        PYTHONPATH="$python_root/src${PYTHONPATH:+:$PYTHONPATH}" \
        "$python_bin" -m taxutils "$@"
}

run_rust() {
    TAXUTILS_GLOBALS="$work/rust-globals" "$tu_bin" --threads 3 "$@"
}

# extract: ordered header parsing and atomic output
run_python extract "$fixtures/valid.fasta" \
    -o "$work/python/extracted.txt" --batch-size 2 > "$work/python/extract.stdout"
run_rust extract "$fixtures/valid.fasta" \
    -o "$work/rust/extracted.txt" --batch-size 2 > "$work/rust/extract.stdout"
cmp "$work/python/extracted.txt" "$work/rust/extracted.txt"

# clean: explicit output and in-place replacement
run_python clean -i "$fixtures/valid.fasta" -o "$work/python/clean.fasta"
run_rust clean -i "$fixtures/valid.fasta" -o "$work/rust/clean.fasta"
cmp "$work/python/clean.fasta" "$work/rust/clean.fasta"

cp "$fixtures/valid.fasta" "$work/python/in-place.fasta"
cp "$fixtures/valid.fasta" "$work/rust/in-place.fasta"
run_python clean -i "$work/python/in-place.fasta"
run_rust clean -i "$work/rust/in-place.fasta"
cmp "$work/python/in-place.fasta" "$work/rust/in-place.fasta"

# grep: versioned, unversioned, missing headers, and byte-bounded batches
run_python grep -i "$fixtures/input.fasta" -a "$fixtures/accessions.txt" \
    -o "$work/python/grep.fasta" --batch-size 17 > "$work/python/grep.stdout"
run_rust grep -i "$fixtures/input.fasta" -a "$fixtures/accessions.txt" \
    -o "$work/rust/grep.fasta" --batch-size 17 > "$work/rust/grep.stdout"
cmp "$work/python/grep.fasta" "$work/rust/grep.fasta"
cmp "$work/python/grep.stdout" "$work/rust/grep.stdout"

run_python grep -i "$fixtures/input.fasta" -a NC_000003 \
    -o "$work/python/grep-no-version.fasta" --no-version
run_rust grep -i "$fixtures/input.fasta" -a NC_000003 \
    -o "$work/rust/grep-no-version.fasta" --no-version
cmp "$work/python/grep-no-version.fasta" "$work/rust/grep-no-version.fasta"

# filter: both modes, SQLite lookup, missing accessions, and record batches
run_python filter -i "$fixtures/input.fasta" -o "$work/python/keep.fasta" \
    --keep-taxids 13 --batch-size 2 > "$work/python/keep.stdout"
run_rust filter -i "$fixtures/input.fasta" -o "$work/rust/keep.fasta" \
    --keep-taxids 13 --batch-size 2 > "$work/rust/keep.stdout"
cmp "$work/python/keep.fasta" "$work/rust/keep.fasta"
cmp "$work/python/keep.stdout" "$work/rust/keep.stdout"

run_python filter -i "$fixtures/input.fasta" -o "$work/python/remove.fasta" \
    --remove-taxids 15 --batch-size 2 > "$work/python/remove.stdout"
run_rust filter -i "$fixtures/input.fasta" -o "$work/rust/remove.fasta" \
    --remove-taxids 15 --batch-size 2 > "$work/rust/remove.stdout"
cmp "$work/python/remove.fasta" "$work/rust/remove.fasta"
cmp "$work/python/remove.stdout" "$work/rust/remove.stdout"

echo "CLI parity passed: extract, clean, grep, and filter"
