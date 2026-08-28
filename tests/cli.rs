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

#[test]
fn seq_selects_an_explicit_model() {
    let input = std::env::temp_dir().join(format!(
        "arpeggia-cli-sequence-model-{}.pdb",
        std::process::id()
    ));
    std::fs::write(
        &input,
        "MODEL        7\n\
         ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
         ENDMDL\n\
         MODEL        9\n\
         ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
         ENDMDL\nEND\n",
    )
    .unwrap();
    let output = arpeggia()
        .args(["seq", "--model", "9", input.to_str().unwrap()])
        .output()
        .unwrap();

    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stdout).contains("A: G"));
}

#[test]
fn sc_calculation_failure_exits_without_a_score() {
    let input = std::env::temp_dir().join(format!(
        "arpeggia-cli-unsupported-radius-{}.pdb",
        std::process::id()
    ));
    std::fs::write(
        &input,
        "ATOM      1  QQ  ALA A   1       0.000   0.000   0.000  1.00 20.00          RN  \n\
         ATOM      2  CB  ALA B   1       3.000   0.000   0.000  1.00 20.00           C  \n\
         END\n",
    )
    .unwrap();
    let output = arpeggia()
        .args(["sc", "--input", input.to_str().unwrap(), "--groups", "A/B"])
        .output()
        .unwrap();

    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("van der Waals radius"));
    assert!(!stderr.contains("INFO arpeggia::cli::sc: SC:"));
}

#[test]
fn rmsd_prints_one_scalar() {
    let input = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
    let output = arpeggia()
        .args(["rmsd", &input, &input, "--residues", "A:1-20"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let value = String::from_utf8(output.stdout)
        .unwrap()
        .trim()
        .parse::<f64>()
        .unwrap();
    assert!(value < 1e-12);
}

#[test]
fn cluster_structs_saves_and_reuses_pairwise_rmsd() {
    let source = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
    let root = std::env::temp_dir().join(format!("arpeggia-cli-clustering-{}", std::process::id()));
    let input = root.join("input");
    let output = root.join("output");
    std::fs::create_dir_all(&input).unwrap();
    for id in ["a", "b", "c"] {
        std::fs::copy(&source, input.join(format!("{id}.pdb"))).unwrap();
    }
    let arguments = [
        "cluster-structs",
        "--input",
        input.to_str().unwrap(),
        "--output",
        output.to_str().unwrap(),
        "--num-clusters",
        "1",
        "--pairwise-rmsd",
        "--residues",
        "A:1-20",
        "--num-threads",
        "2",
    ];
    let first = arpeggia().args(arguments).output().unwrap();
    assert!(
        first.status.success(),
        "{}",
        String::from_utf8_lossy(&first.stderr)
    );
    assert!(output.join("clusters.csv").is_file());
    assert!(output.join("pairwise_rmsd.csv").is_file());

    let second = arpeggia()
        .env("RUST_LOG", "debug")
        .args(arguments)
        .output()
        .unwrap();
    assert!(second.status.success());
    assert!(String::from_utf8_lossy(&second.stderr).contains("Reusing pairwise RMSD cache"));
}

#[test]
fn cluster_structs_rejects_colliding_outputs() {
    let source = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
    let root = std::env::temp_dir().join(format!(
        "arpeggia-cli-cluster-output-collision-{}",
        std::process::id()
    ));
    let input = root.join("input");
    std::fs::create_dir_all(&input).unwrap();
    for id in ["a", "b", "c"] {
        std::fs::copy(&source, input.join(format!("{id}.pdb"))).unwrap();
    }
    let result = arpeggia()
        .args([
            "cluster-structs",
            "--input",
            input.to_str().unwrap(),
            "--output",
            root.to_str().unwrap(),
            "--num-clusters",
            "1",
            "--pairwise-rmsd",
            "--pairwise-filename",
            "clusters",
        ])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(String::from_utf8_lossy(&result.stderr).contains("output paths must differ"));
}

#[test]
fn cluster_structs_preserves_pairwise_work_and_rejects_bad_cache() {
    let source = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
    let root = std::env::temp_dir().join(format!(
        "arpeggia-cli-cluster-failure-{}",
        std::process::id()
    ));
    let input = root.join("input");
    let output = root.join("output");
    std::fs::create_dir_all(&input).unwrap();
    for id in ["a", "b", "c"] {
        std::fs::copy(&source, input.join(format!("{id}.pdb"))).unwrap();
    }
    let failed = arpeggia()
        .args([
            "cluster-structs",
            "--input",
            input.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--num-clusters",
            "1",
            "--max-iterations",
            "1",
            "--pairwise-rmsd",
            "--residues",
            "A:1-20",
        ])
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert!(output.join("pairwise_rmsd.csv").is_file());
    assert!(!output.join("clusters.csv").is_file());

    std::fs::write(output.join("pairwise_rmsd.csv"), "broken\ncache\n").unwrap();
    let malformed = arpeggia()
        .args([
            "cluster-structs",
            "--input",
            input.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--num-clusters",
            "1",
            "--pairwise-rmsd",
        ])
        .output()
        .unwrap();
    assert!(!malformed.status.success());
    assert!(String::from_utf8_lossy(&malformed.stderr).contains("remove it to recalculate"));
}
