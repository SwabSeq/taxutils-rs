use std::path::Path;
use std::process::Command;

#[test]
fn cli_outputs_match_python_taxutils() {
    let manifest = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("bash")
        .arg(manifest.join("tests/cli_parity.sh"))
        .env("TU_BIN", env!("CARGO_BIN_EXE_tu"))
        .output()
        .expect("run CLI parity script");
    assert!(
        output.status.success(),
        "CLI parity failed.\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}
