use std::collections::HashSet;
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use taxutils::{CancellationToken, FilterMode};

const FASTA: &[u8] = b">NC_000001.1 first\r\nAC\r\nGT\r\n>NC_000002.1 other\nTT\n>NC_000001.1 duplicate\nCC\n>NC_000003.2 last\nGG";

fn encode(data: &[u8], format: &str) -> Vec<u8> {
    match format {
        "gz" => {
            let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
            encoder.write_all(data).unwrap();
            encoder.finish().unwrap()
        }
        "zst" => {
            let mut encoder = zstd::stream::write::Encoder::new(Vec::new(), 3).unwrap();
            encoder.include_checksum(true).unwrap();
            encoder.write_all(data).unwrap();
            encoder.finish().unwrap()
        }
        _ => data.to_vec(),
    }
}

fn decode(path: &Path, format: &str) -> Vec<u8> {
    let data = fs::read(path).unwrap();
    match format {
        "gz" => {
            assert!(data.starts_with(&[0x1f, 0x8b]));
            let mut output = Vec::new();
            MultiGzDecoder::new(&data[..])
                .read_to_end(&mut output)
                .unwrap();
            output
        }
        "zst" => {
            assert!(data.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]));
            zstd::stream::decode_all(&data[..]).unwrap()
        }
        _ => data,
    }
}

fn globals(directory: &Path) -> PathBuf {
    let target = directory.join("globals");
    fs::create_dir(&target).unwrap();
    let fixtures = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/ncbi");
    for name in ["names.dmp", "nodes.dmp", "targets.json"] {
        fs::copy(fixtures.join(name), target.join(name)).unwrap();
    }
    for source in ["nucl_gb", "nucl_wgs"] {
        fs::write(
            target.join(format!("{source}.accession2taxid.gz")),
            encode(
                &fs::read(fixtures.join(format!("{source}.accession2taxid.tsv"))).unwrap(),
                "gz",
            ),
        )
        .unwrap();
    }
    target
}

fn run(
    command: &str,
    input: &Path,
    output: Option<&Path>,
    cache: &Path,
    query: Option<&Path>,
) -> Output {
    let mut cli = Command::new(env!("CARGO_BIN_EXE_tu"));
    cli.env("TAXUTILS_GLOBALS", cache)
        .args(["--threads", "2", command]);
    if command != "extract" {
        cli.arg("-i");
    }
    cli.arg(input);
    if let Some(output) = output {
        cli.arg("-o").arg(output);
    }
    if command == "grep" {
        cli.arg("-a");
        if let Some(query) = query {
            cli.arg(query);
        } else {
            cli.arg("NC_000001.1");
        }
    }
    if command == "filter" {
        cli.arg("--keep-taxids");
        if let Some(query) = query {
            cli.arg(query);
        } else {
            cli.arg("13");
        }
    }
    if ["extract", "filter", "grep"].contains(&command) {
        cli.args(["--batch-size", "1"]);
    }
    cli.output().unwrap()
}

fn success(output: &Output) {
    assert!(
        output.status.success(),
        "stdout: {}\nstderr: {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn every_cli_command_accepts_and_writes_each_format() {
    let dir = tempfile::tempdir().unwrap();
    let cache = globals(dir.path());
    let input = dir.path().join("input.fna");
    let output = dir.path().join("output");
    for command in ["extract", "clean", "deduplicate", "grep", "filter"] {
        fs::write(&input, FASTA).unwrap();
        let baseline = run(command, &input, Some(&output), &cache, None);
        success(&baseline);
        let expected = fs::read(&output).unwrap();
        for input_format in ["plain", "gz", "zst"] {
            // Content detection works even with a misleading plain suffix.
            fs::write(&input, encode(FASTA, input_format)).unwrap();
            for (suffix, output_format) in [
                ("txt", "plain"),
                ("GZ", "gz"),
                ("ZST", "zst"),
                ("ZsTd", "zst"),
            ] {
                let destination = dir.path().join(format!("output.{suffix}"));
                let result = run(command, &input, Some(&destination), &cache, None);
                success(&result);
                assert_eq!(
                    decode(&destination, output_format),
                    expected,
                    "{command}: {input_format} -> {suffix}"
                );
                if command != "extract" {
                    assert_eq!(result.stdout, baseline.stdout);
                }
            }
        }
    }
}

#[test]
fn rewrites_preserve_content_format_and_explicit_paths_use_suffix() {
    let dir = tempfile::tempdir().unwrap();
    let cache = globals(dir.path());
    let baseline_input = dir.path().join("baseline.fna");
    let baseline_output = dir.path().join("baseline-output.fna");
    for command in ["clean", "deduplicate", "filter"] {
        fs::write(&baseline_input, FASTA).unwrap();
        success(&run(
            command,
            &baseline_input,
            Some(&baseline_output),
            &cache,
            None,
        ));
        let expected = fs::read(&baseline_output).unwrap();
        for format in ["gz", "zst"] {
            for suffix in ["fna", format] {
                let input = dir.path().join(format!("rewrite.{suffix}"));
                fs::write(&input, encode(FASTA, format)).unwrap();
                success(&run(command, &input, None, &cache, None));
                assert_eq!(decode(&input, format), expected);
                fs::write(&input, encode(FASTA, format)).unwrap();
                success(&run(command, &input, Some(&input), &cache, None));
                assert_eq!(
                    decode(&input, if suffix == "fna" { "plain" } else { format }),
                    expected
                );
            }
        }
    }
    // grep must open input before installing same-path output.
    let input = dir.path().join("grep.zst");
    fs::write(&input, encode(FASTA, "gz")).unwrap();
    success(&run("grep", &input, Some(&input), &cache, None));
    assert_eq!(
        decode(&input, "zst"),
        b">NC_000001.1 first\r\nAC\r\nGT\r\n>NC_000001.1 duplicate\nCC\n"
    );
}

#[test]
fn concatenated_streams_skippable_frames_and_query_files() {
    let dir = tempfile::tempdir().unwrap();
    let cache = globals(dir.path());
    let input = dir.path().join("input");
    let output = dir.path().join("output");
    let query = dir.path().join("query.txt");
    // Split in the middle of a line to ensure decoder boundaries are invisible.
    let split = 20;
    for format in ["gz", "zst"] {
        let mut stream = encode(&FASTA[..split], format);
        if format == "zst" {
            // Skippable-frame magic, payload size, then payload.
            let skip = [0x50, 0x2a, 0x4d, 0x18, 3, 0, 0, 0, 1, 2, 3];
            stream.extend(skip);
            let mut prefixed = skip.to_vec();
            prefixed.extend(stream);
            stream = prefixed;
        }
        stream.extend(encode(&FASTA[split..], format));
        fs::write(&input, stream).unwrap();
        success(&run("extract", &input, Some(&output), &cache, None));
        assert_eq!(
            fs::read(&output).unwrap(),
            b"NC_000001.1\nNC_000002.1\nNC_000001.1\nNC_000003.2\n"
        );
        // Feed compressed extract output straight back into grep.
        let query_output = dir.path().join(format!("accessions.{format}"));
        success(&run("extract", &input, Some(&query_output), &cache, None));
        success(&run(
            "grep",
            &input,
            Some(&output),
            &cache,
            Some(&query_output),
        ));
        assert_eq!(fs::read(&output).unwrap(), FASTA);
        fs::write(&query, encode(b"13\n", format)).unwrap();
        success(&run("filter", &input, Some(&output), &cache, Some(&query)));
        assert!(
            fs::read(&output)
                .unwrap()
                .starts_with(b">NC_000001.1 first")
        );
    }
}

#[test]
fn corrupt_and_truncated_inputs_preserve_destinations_and_remove_temporary_files() {
    let dir = tempfile::tempdir().unwrap();
    let cache = globals(dir.path());
    let input = dir.path().join("input");
    let output = dir.path().join("output.zst");
    for format in ["gz", "zst"] {
        let bytes = encode(FASTA, format);
        let mut corrupted = bytes.clone();
        let last = corrupted.len() - 1;
        corrupted[last] ^= 0xff;
        for invalid in [bytes[..bytes.len() - 1].to_vec(), corrupted] {
            for command in ["extract", "clean", "deduplicate", "grep", "filter"] {
                fs::write(&input, &invalid).unwrap();
                fs::write(&output, b"existing destination").unwrap();
                let result = run(command, &input, Some(&output), &cache, None);
                assert!(
                    !result.status.success(),
                    "{command} accepted malformed {format}"
                );
                assert!(
                    String::from_utf8_lossy(&result.stderr).contains(&input.display().to_string())
                );
                assert_eq!(fs::read(&output).unwrap(), b"existing destination");
                assert_eq!(fs::read(&input).unwrap(), invalid);
                assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 3);
            }
            for command in ["clean", "deduplicate", "filter"] {
                let result = run(command, &input, None, &cache, None);
                assert!(!result.status.success());
                assert_eq!(fs::read(&input).unwrap(), invalid);
                assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 3);
            }
        }
    }
}

#[test]
fn empty_compressed_input_is_valid_for_all_commands() {
    let dir = tempfile::tempdir().unwrap();
    let cache = globals(dir.path());
    let input = dir.path().join("input");
    let output = dir.path().join("output.gz");
    for format in ["plain", "gz", "zst"] {
        fs::write(&input, encode(b"", format)).unwrap();
        for command in ["extract", "clean", "deduplicate", "grep", "filter"] {
            success(&run(command, &input, Some(&output), &cache, None));
            assert!(decode(&output, "gz").is_empty());
        }
    }
}

#[test]
fn cancelled_library_operations_preserve_compressed_input_and_destination() {
    let dir = tempfile::tempdir().unwrap();
    let cache = globals(dir.path());
    let input = dir.path().join("input.zst");
    let output = dir.path().join("output.gz");
    let original = encode(FASTA, "zst");
    fs::write(&input, &original).unwrap();
    fs::write(&output, b"existing").unwrap();
    let cancellation = CancellationToken::default();
    cancellation.cancel();
    assert!(
        taxutils::extract_accessions_with_cancel(&input, &output, 1, Some(1), &cancellation)
            .is_err()
    );
    assert!(
        taxutils::grep_fasta_with_cancel(
            &input,
            "NC_000001.1",
            &output,
            true,
            1,
            false,
            Some(1),
            &cancellation
        )
        .is_err()
    );
    for destination in [None, Some(output.as_path())] {
        assert!(
            taxutils::clean_fasta_headers_with_cancel(
                &input,
                destination,
                false,
                Some(1),
                &cancellation
            )
            .is_err()
        );
        assert!(
            taxutils::deduplicate_fasta_with_cancel(&input, destination, Some(1), &cancellation)
                .is_err()
        );
        assert!(
            taxutils::filter_fasta_with_options_and_cancel(
                &input,
                destination,
                &HashSet::from([13]),
                FilterMode::Keep,
                1,
                false,
                &cache,
                false,
                Some(1),
                &cancellation
            )
            .is_err()
        );
    }
    assert_eq!(fs::read(&input).unwrap(), original);
    assert_eq!(fs::read(&output).unwrap(), b"existing");
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 3);
}

#[cfg(unix)]
#[test]
fn atomic_rewrites_retain_destination_permissions() {
    use std::os::unix::fs::PermissionsExt;
    let dir = tempfile::tempdir().unwrap();
    let cache = globals(dir.path());
    let input = dir.path().join("input.zst");
    let output = dir.path().join("output.gz");
    for command in ["extract", "clean", "deduplicate", "grep", "filter"] {
        fs::write(&input, encode(FASTA, "zst")).unwrap();
        fs::write(&output, b"old").unwrap();
        fs::set_permissions(&output, fs::Permissions::from_mode(0o640)).unwrap();
        success(&run(command, &input, Some(&output), &cache, None));
        assert_eq!(
            fs::metadata(&output).unwrap().permissions().mode() & 0o777,
            0o640
        );
    }
}
