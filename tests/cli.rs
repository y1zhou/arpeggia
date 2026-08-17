use std::process::Command;

fn arpeggia() -> Command {
    Command::new(env!("CARGO_BIN_EXE_arpeggia"))
}

#[test]
fn missing_input_exits_nonzero() {
    let output = arpeggia()
        .args(["seq", "this-file-does-not-exist.pdb"])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("I/O error"));
}

#[test]
fn unsafe_output_filename_exits_nonzero() {
    let input = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
    let output_dir =
        std::env::temp_dir().join(format!("arpeggia-cli-output-{}", std::process::id()));
    let output = arpeggia()
        .args([
            "sasa",
            "--input",
            &input,
            "--output",
            output_dir.to_str().unwrap(),
            "--filename",
            "../escape",
            "--num-points",
            "1",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("one normal path component"));
}

#[test]
fn nonfinite_scientific_parameter_exits_nonzero() {
    let input = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
    let output_dir =
        std::env::temp_dir().join(format!("arpeggia-cli-parameter-{}", std::process::id()));
    let output = arpeggia()
        .args([
            "sasa",
            "--input",
            &input,
            "--output",
            output_dir.to_str().unwrap(),
            "--probe-radius",
            "NaN",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("finite and positive"));
}

#[test]
fn ringless_contacts_succeed_through_the_cli() {
    let input =
        std::env::temp_dir().join(format!("arpeggia-cli-ringless-{}.pdb", std::process::id()));
    let output_dir =
        std::env::temp_dir().join(format!("arpeggia-cli-ringless-{}", std::process::id()));
    std::fs::write(
        &input,
        "ATOM      1  NZ  LYS A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
         ATOM      2  OD1 ASP B   1       2.500   0.000   0.000  1.00 20.00           O  \n\
         END\n",
    )
    .unwrap();
    let output = arpeggia()
        .args([
            "contacts",
            "--input",
            input.to_str().unwrap(),
            "--output",
            output_dir.to_str().unwrap(),
            "--groups",
            "A/B",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn conformer_selection_is_warned_on_stderr() {
    let input =
        std::env::temp_dir().join(format!("arpeggia-cli-conformer-{}.pdb", std::process::id()));
    std::fs::write(
        &input,
        "ATOM      1  CB AALA A   1       0.000   0.000   0.000  0.50 20.00           C  \n\
         ATOM      2  CB BALA A   1       1.000   0.000   0.000  0.50 20.00           C  \n\
         END\n",
    )
    .unwrap();
    let output = arpeggia()
        .args(["seq", input.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("CONFORMER_SELECTED"));
}
