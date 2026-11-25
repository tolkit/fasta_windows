//! End-to-end regression tests: run the binary on a small FASTA and compare
//! every output file byte for byte with tests/data/expected/<variant>/.
//!
//! The expected files were produced by fasta_windows 0.2.6 (the HashMap-based
//! implementation), so these tests pin the optimised code to the original
//! output. The fixture has no IUPAC codes, where 0.2.6 produced ragged k-mer
//! rows, and its record IDs are already sorted, so 0.2.6's sort-by-ID order
//! equals FASTA order.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn data_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
}

/// Run fasta_windows in a fresh temporary directory; return its fw_out/.
fn run(name: &str, input: &str, args: &[&str]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("fw_golden_{}_{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_fasta_windows"))
        .current_dir(&dir)
        .arg("-f")
        .arg(data_dir().join(input))
        .args(["-o", "x"])
        .args(args)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "fasta_windows failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    dir.join("fw_out")
}

fn assert_matches_expected(variant: &str, out: &Path) {
    let expected_dir = data_dir().join("expected").join(variant);
    let mut checked = 0;
    for entry in fs::read_dir(&expected_dir).unwrap() {
        let entry = entry.unwrap();
        let name = entry.file_name();
        let expected = fs::read(entry.path()).unwrap();
        let got = fs::read(out.join(&name))
            .unwrap_or_else(|e| panic!("missing output {}: {e}", name.to_string_lossy()));
        if expected != got {
            let first_diff = String::from_utf8_lossy(&expected)
                .lines()
                .zip(String::from_utf8_lossy(&got).lines())
                .enumerate()
                .find(|(_, (e, g))| e != g)
                .map(|(i, (e, g))| format!("line {}:\nexpected: {e}\n     got: {g}", i + 1))
                .unwrap_or_else(|| "outputs differ in length".to_string());
            panic!(
                "{variant}/{} differs. First difference at {first_diff}",
                name.to_string_lossy()
            );
        }
        checked += 1;
    }
    assert!(checked > 0, "no expected files for {variant}");
    assert_eq!(
        fs::read_dir(out).unwrap().count(),
        checked,
        "unexpected extra output files"
    );
    let _ = fs::remove_dir_all(out.parent().unwrap());
}

#[test]
fn golden_default() {
    let out = run("default", "small.fa", &["-w", "1000"]);
    assert_matches_expected("default", &out);
}

#[test]
fn golden_ctw_description() {
    let out = run("ctw_description", "small.fa", &["-c", "-d", "-w", "333"]);
    assert_matches_expected("ctw_description", &out);
}

#[test]
fn golden_masked() {
    let out = run("masked", "small.fa", &["-m", "-w", "500"]);
    assert_matches_expected("masked", &out);
}

#[test]
fn golden_gzip_input() {
    let out = run("gzip", "small.fa.gz", &["-c", "-d", "-w", "333"]);
    assert_matches_expected("ctw_description", &out);
}

#[cfg(feature = "entropy")]
#[test]
fn golden_entropy() {
    let out = run("entropy", "small.fa", &["-e", "-w", "400"]);
    assert_matches_expected("entropy", &out);
}
