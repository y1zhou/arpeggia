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
fn rmsd_rejects_selection_before_structure_io() {
    let output = arpeggia()
        .args([
            "rmsd",
            "missing-reference.pdb",
            "missing-query.pdb",
            "--superpose-residues",
            "A:",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(stderr.contains("Superposition Selection"));
    assert!(stderr.contains("empty residue selection"));

    let output = arpeggia()
        .args([
            "rmsd",
            "missing-reference.pdb",
            "missing-query.pdb",
            "--rmsd-residues",
            "A:",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("RMSD Selection"));
}

#[test]
fn cluster_structs_rejects_options_before_directory_io() {
    let output = arpeggia()
        .args([
            "cluster-structs",
            "--input",
            "missing-structure-directory",
            "--output",
            "unused-output",
        ])
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(
        String::from_utf8_lossy(&output.stderr)
            .contains("one of num_clusters or max_clusters is required")
    );
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
fn rmsd_reports_detailed_json() {
    let input = format!("{}/test-data/1ubq.pdb", env!("CARGO_MANIFEST_DIR"));
    let output = arpeggia()
        .args([
            "rmsd",
            &input,
            &input,
            "--superpose-residues",
            "A:1-20",
            "--rmsd-residues",
            "A:1-20",
            "--json",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(result["rmsd"].as_f64().unwrap() < 1e-12);
    assert_eq!(result["core_rmsd"], 0.0);
    assert_eq!(result["evaluation_atoms"], 20);
    assert_eq!(result["retained_fit_atoms"], 20);
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
        "--superpose-residues",
        "A:1-20",
        "--rmsd-residues",
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
        .args([
            "cluster-structs",
            "--input",
            input.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--num-clusters",
            "1",
            "--pairwise-rmsd",
            "--superpose-residues",
            "A:1-10",
            "--rmsd-residues",
            "A:1-10",
        ])
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
    std::fs::create_dir_all(output.join("clusters.csv")).unwrap();
    let failed = arpeggia()
        .args([
            "cluster-structs",
            "--input",
            input.to_str().unwrap(),
            "--output",
            output.to_str().unwrap(),
            "--num-clusters",
            "1",
            "--pairwise-rmsd",
            "--superpose-residues",
            "A:1-20",
            "--rmsd-residues",
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

#[test]
fn sequence_cli_reports_metrics_and_empty_local_json() {
    let result = arpeggia()
        .args([
            "align-seqs",
            "GGACDEFGHIKGG",
            "ACDEFGHIK",
            "--mode",
            "semi-global",
            "--json",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let data: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(data["reference_span"], serde_json::json!([2, 11]));
    assert_eq!(data["identity_shorter"], 1.0);
    assert_eq!(data["edit_distance"], 4);
    let result = arpeggia()
        .args(["align-seqs", "AAAA", "WWWW", "--mode", "local", "--json"])
        .output()
        .unwrap();
    assert!(result.status.success());
    let data: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(data["identity_alignment"], serde_json::Value::Null);
    assert_eq!(data["coverage_shorter"], 0.0);
}

#[test]
fn alignment_display_flags_preserve_plain_json_and_wrapping() {
    let run = |extra: &[&str]| {
        arpeggia()
            .args(["align-seqs", "ACDEFGHIKLMN", "ACDYGHIKLMN"])
            .args(extra)
            .output()
            .unwrap()
    };
    let plain = run(&["--width", "40", "--no-rulers"]);
    assert!(plain.status.success());
    let text = String::from_utf8(plain.stdout).unwrap();
    assert!(!text.contains('\x1b'));
    assert!(!text.contains("operations"));
    assert!(text.lines().all(|line| line.len() <= 40));
    let ruled = run(&["--width", "40"]);
    assert_eq!(
        String::from_utf8(ruled.stdout).unwrap().lines().count(),
        text.lines().count() + 2
    );
    let colored = run(&["--color", "always"]);
    let text = String::from_utf8(colored.stdout).unwrap();
    assert!(text.contains("\x1b[31m") && text.contains("\x1b[34m:\x1b[0m"));
    let json = run(&["--json", "--color", "always"]);
    let data: serde_json::Value = serde_json::from_slice(&json.stdout).unwrap();
    assert_eq!(data["reference_name"], "Reference");
    assert_eq!(data["query_name"], "Query");
    let named = run(&["--reference-name", "Wild type", "--query-name", "Mutant"]);
    assert!(named.status.success());
    let text = String::from_utf8(named.stdout).unwrap();
    assert!(text.contains("Wild type") && text.contains("Mutant"));
    let named_json = run(&[
        "--reference-name",
        "Wild type",
        "--query-name",
        "Mutant",
        "--json",
    ]);
    let named_data: serde_json::Value = serde_json::from_slice(&named_json.stdout).unwrap();
    assert_eq!(named_data["reference_name"], "Wild type");
    assert_eq!(named_data["query_name"], "Mutant");
    assert_eq!(named_data["score"], data["score"]);
    assert_eq!(data["aligned_query"], "ACD-YGHIKLMN");
    assert_eq!(data["operations"], "   -:       ");
    assert_eq!(data["mismatches"], 1);
    assert!(data.get("columns").is_none());
    assert!(!run(&["--width", "1"]).status.success());
}

const ANTIBODY_SEQUENCE: &str = "QVQLVQSGAEVKRPGSSVTVSCKASGGSFSTYALSWVRQAPGRGLEWMGGVIPLLTITNYAPRFQGRITITADRSTSTAYLELNSLRPEDTAVYYCAREGTTGKPIGAFAHWGQGTLVTVSS";

#[test]
fn antibody_cli_preserves_names_and_structured_numbering() {
    let output = arpeggia()
        .args([
            "number-antibody",
            ANTIBODY_SEQUENCE,
            "--name",
            "WT",
            "--scheme",
            "chothia",
            "--species",
            "rat,rabbit",
            "--json",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(result["name"], "WT");
    assert_eq!(result["scheme"], "martin");
    assert_eq!(result["cdr_definition"], "martin");
    assert!(result["residues"].as_array().unwrap().len() > 100);
    for segment in ["v_match", "j_match"] {
        for hit in result[segment]["hits"].as_array().unwrap() {
            for reference in hit["references"].as_array().unwrap() {
                let species = reference["species"].as_str().unwrap();
                assert!(
                    species.starts_with("Rattus norvegicus")
                        || species.starts_with("Oryctolagus cuniculus")
                );
            }
        }
    }
    assert!(!output.stdout.contains(&0x1b));
    let output = arpeggia()
        .args([
            "align-antibodies",
            ANTIBODY_SEQUENCE,
            ANTIBODY_SEQUENCE,
            ANTIBODY_SEQUENCE,
            "--names",
            "WT,,Mutant",
            "--reference-index",
            "2",
            "--json",
        ])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(result["reference_index"], 2);
    assert_eq!(result["antibodies"][0]["name"], "WT");
    assert_eq!(result["antibodies"][1]["name"], "Seq002");
    assert_eq!(result["antibodies"][2]["name"], "Mutant");
}

#[test]
fn antibody_cli_ruler_flag_preserves_compact_reference_first_blocks() {
    for command in ["number-antibody", "align-antibodies"] {
        let run = |rulers: bool| {
            let mut args = vec![command, ANTIBODY_SEQUENCE];
            if command == "number-antibody" {
                args.extend(["--name", "WT"]);
            } else {
                args.extend([
                    ANTIBODY_SEQUENCE,
                    "--names",
                    "Mutant,WT",
                    "--reference-index",
                    "1",
                ]);
            }
            args.extend(["--width", "80", "--color", "never"]);
            if !rulers {
                args.push("--no-rulers");
            }
            let output = arpeggia().args(args).output().unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            String::from_utf8(output.stdout).unwrap()
        };
        let full = run(true);
        let compact = run(false);
        assert!(full.contains("Reference: WT"));
        assert!(!full.contains('\x1b'));
        assert!(full.lines().all(|l| l.len() <= 80));
        let full_blocks: Vec<_> = full.split("\n\n").skip(1).collect();
        let compact_blocks: Vec<_> = compact.split("\n\n").skip(1).collect();
        assert!(!full_blocks.is_empty());
        assert_eq!(full_blocks.len(), compact_blocks.len());
        for (full, compact) in full_blocks.iter().zip(compact_blocks) {
            let rows: Vec<_> = full.lines().collect();
            assert_eq!(rows.len(), 6);
            assert!(rows[2].starts_with("WT "));
            assert_eq!(
                compact.lines().collect::<Vec<_>>(),
                [rows[0], rows[2], rows[4], rows[5]]
            );
        }
    }
}

#[test]
fn antibody_cli_imputation_summary_requires_impute() {
    for command in ["number-antibody", "align-antibodies"] {
        for impute in [false, true] {
            let mut args = vec![command, ANTIBODY_SEQUENCE, "--color", "never"];
            if impute {
                args.push("--impute");
            }
            let output = arpeggia().args(args).output().unwrap();
            assert!(output.status.success(), "{:?}", output);
            let text = String::from_utf8(output.stdout).unwrap();
            assert_eq!(text.contains("imputed residues:"), impute);
            if impute {
                assert!(text.contains("imputed residues: 0"));
            }
            assert!(text.contains("CDR regions:"));
            assert!(!text.contains("Total CDR"));
        }
    }
}

#[test]
fn antibody_cli_can_skip_germlines_without_disabling_numbering() {
    for command in ["number-antibody", "align-antibodies"] {
        let args = [
            command,
            ANTIBODY_SEQUENCE,
            "--no-germlines",
            "--species",
            "rat,rabbit",
        ];
        let output = arpeggia().args(args).arg("--json").output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty());
        let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        let result = if command == "number-antibody" {
            &value
        } else {
            &value["antibodies"][0]
        };
        assert_eq!(result["germlines_searched"], false);
        assert!(result["v_match"].is_null() && result["j_match"].is_null());
        assert!(result["diagnostics"].as_array().unwrap().is_empty());
        assert_eq!(
            result["residues"].as_array().unwrap().len(),
            ANTIBODY_SEQUENCE.len()
        );
        let conflict = arpeggia().args(args).arg("--impute").output().unwrap();
        assert!(!conflict.status.success());
        let error = String::from_utf8_lossy(&conflict.stderr);
        assert!(error.contains("--no-germlines") && error.contains("--impute"));
    }
}

#[test]
fn antibody_cli_escapes_names_in_diagnostics() {
    let output = arpeggia()
        .args([
            "align-antibodies",
            &ANTIBODY_SEQUENCE[5..],
            "--names",
            "partial\nforged\rname",
            "--no-germlines",
            "--color",
            "never",
        ])
        .output()
        .unwrap();
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(output.status.success(), "{stderr}");
    assert!(
        stderr.contains(r#""partial\nforged\rname": PARTIAL_DOMAIN"#),
        "{stderr}"
    );
    assert_eq!(stderr.lines().count(), 1, "{stderr}");
    assert!(!stderr.contains('\r'), "{stderr}");
}

#[test]
fn antibody_cli_rejects_invalid_names_and_conventions() {
    for options in [
        vec!["align-antibodies", ANTIBODY_SEQUENCE, "--names", "one,two"],
        vec![
            "number-antibody",
            ANTIBODY_SEQUENCE,
            "--cdr-definition",
            "chothia",
        ],
        vec![
            "align-antibodies",
            ANTIBODY_SEQUENCE,
            "--reference-index",
            "1",
        ],
    ] {
        let output = arpeggia().args(options).output().unwrap();
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("invalid argument"));
    }
}
