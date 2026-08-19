//! Contact analysis functionality for protein structures.
//!
//! This module provides functions for analyzing atomic and ring contacts
//! in PDB and mmCIF structures.
pub mod aromatic;
pub mod complex;
pub mod hbond;
pub mod hydrophobic;
pub mod ionic;
pub mod residues;
pub mod structs;
pub mod vdw;

// Re-exports
pub use aromatic::{find_cation_pi, find_pi_pi};
pub use complex::*;
pub use hbond::{find_hydrogen_bond, find_weak_hydrogen_bond};
pub use hydrophobic::find_hydrophobic_contact;
pub use ionic::{find_ionic_bond_with_protonation, find_ionic_repulsion_with_protonation};
pub use residues::ResidueId;
pub use structs::*;
pub use vdw::find_vdw_contact;

use pdbtbx::*;
use polars::prelude::*;
use std::collections::HashMap;
use tracing::{debug, warn};

/// Analyze contacts with explicit connectivity and histidine-protonation policy.
#[allow(clippy::too_many_arguments)]
pub fn analyze_contacts(
    pdb: &PDB,
    metadata: Option<&crate::StructureMetadata>,
    groups: &str,
    vdw_comp: f64,
    dist_cutoff: f64,
    protonation: ProtonationMode,
    ph: f64,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    if !vdw_comp.is_finite() || vdw_comp < 0.0 {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "vdw_comp must be finite and non-negative".into(),
        ));
    }
    if !dist_cutoff.is_finite() || dist_cutoff <= 0.0 {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "dist_cutoff must be finite and positive".into(),
        ));
    }
    if !ph.is_finite() {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "pH must be finite".into(),
        ));
    }
    let mut prepared = pdb.clone();
    let mut warnings = crate::structure::select_conformers(&mut prepared);
    let (complex, ring_warnings) = InteractionComplex::new_with_options(
        &prepared,
        metadata,
        groups,
        vdw_comp,
        dist_cutoff,
        protonation,
        ph,
    )?;
    warnings.extend(ring_warnings.into_iter().map(|message| {
        crate::AnalysisWarning::new(crate::WarningCode::IncompleteGeometry, message)
    }));
    warnings.extend(protonation_warnings(&prepared, protonation));
    let missing_donors = hbond::count_donors_without_explicit_hydrogen(&prepared, metadata);
    if missing_donors > 0 {
        warnings.push(crate::AnalysisWarning::new(
            crate::WarningCode::MissingDonorHydrogen,
            format!(
                "{missing_donors} donor atoms lack an associated explicit hydrogen; \
                 distance-only polar contacts remain available"
            ),
        ));
    }
    Ok(crate::Analysis::new(
        contacts_from_complex(complex, Vec::new()),
        warnings,
    ))
}

fn protonation_warnings(pdb: &PDB, mode: ProtonationMode) -> Vec<crate::AnalysisWarning> {
    pdb.residues()
        .filter(|residue| ionic::is_histidine(residue.name().unwrap_or("")))
        .filter(|residue| ionic::explicit_histidine_charge(residue).is_none())
        .map(|residue| {
            crate::AnalysisWarning::new(
                crate::WarningCode::UnresolvedHistidine,
                format!(
                    "histidine {} has no explicit protonation evidence; policy {mode:?} applied",
                    residue.serial_number()
                ),
            )
        })
        .collect()
}

/// Calculate atomic and ring contacts in a PDB structure.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `groups` - Chain groups specification (e.g., "A,B/C,D" or "/" for all-to-all)
/// * `vdw_comp` - VdW radii compensation factor (typically 0.1)
/// * `dist_cutoff` - Distance cutoff for neighbor searches (typically 6.5 Å)
///
/// # Returns
///
/// A Polars `DataFrame` containing all identified contacts with columns:
/// - model, interaction, distance
/// - `from_chain`, `from_resn`, `from_resi`, `from_insertion`, `from_alt_loc`, `from_atomn`, `from_atomi`
/// - `to_chain`, `to_resn`, `to_resi`, `to_insertion`, `to_altloc`, `to_atomn`, `to_atomi`
/// - `sc_centroid_dist`, `sc_dihedral`, `sc_centroid_angle`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_contacts};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let pdb = load_model(&input_file).unwrap().value;
/// let contacts_df = get_contacts(
///     &pdb, "/", 0.1, 6.5,
///     arpeggia::ProtonationMode::AllCharged, 7.4,
/// ).unwrap().value;
/// println!("Found {} contacts", contacts_df.height());
/// ```
pub fn get_contacts(
    pdb: &PDB,
    groups: &str,
    vdw_comp: f64,
    dist_cutoff: f64,
    protonation: ProtonationMode,
    ph: f64,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    analyze_contacts(pdb, None, groups, vdw_comp, dist_cutoff, protonation, ph)
}

/// Calculate contacts using explicit connectivity parsed from the input file.
pub fn get_contacts_with_metadata(
    pdb: &PDB,
    metadata: &crate::StructureMetadata,
    groups: &str,
    vdw_comp: f64,
    dist_cutoff: f64,
    protonation: ProtonationMode,
    ph: f64,
) -> crate::ArpeggiaResult<crate::Analysis<DataFrame>> {
    analyze_contacts(
        pdb,
        Some(metadata),
        groups,
        vdw_comp,
        dist_cutoff,
        protonation,
        ph,
    )
}

fn contacts_from_complex(
    i_complex: InteractionComplex<'_>,
    build_ring_err: Vec<String>,
) -> DataFrame {
    if !build_ring_err.is_empty() {
        for e in &build_ring_err {
            warn!("{e}");
        }
    }

    // Information on the sequence of the chains in the model
    debug!(
        "Parsed ligand chains {lig:?}; receptor chains {receptor:?}",
        lig = i_complex.ligand,
        receptor = i_complex.receptor
    );

    // Find interactions
    let atomic_contacts = i_complex.get_atomic_contacts();
    let df_atomic = results_to_df(&atomic_contacts);
    debug!(
        "Found {} atom-atom contacts\n{}",
        df_atomic.height(),
        df_atomic
    );

    let mut ring_contacts: Vec<ResultEntry> = Vec::new();
    ring_contacts.extend(i_complex.get_ring_atom_contacts());
    ring_contacts.extend(i_complex.get_ring_ring_contacts());
    let df_ring = results_to_df(&ring_contacts);
    debug!("Found {} ring contacts\n{}", df_ring.height(), df_ring);

    // Annotate sidechain centroid distances and dihedrals
    let mut sc_dist_dihedrals = HashMap::new();
    sc_dist_dihedrals.extend(i_complex.collect_sc_stats(&atomic_contacts));
    sc_dist_dihedrals.extend(i_complex.collect_sc_stats(&ring_contacts));

    let df_sc_stats = sc_results_to_df(&sc_dist_dihedrals);

    // Concatenate dataframes
    let contacts_sc_join_cols = [
        "model",
        "from_chain",
        "from_resi",
        "from_insertion",
        "from_altloc",
        "to_chain",
        "to_resi",
        "to_insertion",
        "to_altloc",
    ];
    df_atomic
        .vstack(&df_ring)
        .unwrap()
        .join(
            &df_sc_stats,
            contacts_sc_join_cols,
            contacts_sc_join_cols,
            JoinArgs::new(JoinType::Left),
            None,
        )
        .unwrap()
        .sort(
            [
                "model",
                "from_chain",
                "to_chain",
                "from_resi",
                "from_altloc",
                "from_atomi",
                "to_resi",
                "to_altloc",
                "to_atomi",
                "interaction",
            ],
            SortMultipleOptions::default(),
        )
        .unwrap()
}

/// Convert a slice of `ResultEntry` into a Polars `DataFrame`.
pub(crate) fn results_to_df(res: &[ResultEntry]) -> DataFrame {
    df!(
        "model" => res.iter().map(|x| x.model as u64).collect::<Vec<u64>>(),
        "interaction" => res.iter().map(|x| x.interaction.to_string()).collect::<Vec<String>>(),
        "distance" => res.iter().map(|x| x.distance as f32).collect::<Vec<f32>>(),
        "from_chain" => res.iter().map(|x| x.ligand.chain.clone()).collect::<Vec<String>>(),
        "from_resn" => res.iter().map(|x| x.ligand.resn.clone()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|x| x.ligand.resi as i64).collect::<Vec<i64>>(),
        "from_insertion" => res.iter().map(|x| x.ligand.insertion.clone()).collect::<Vec<String>>(),
        "from_altloc" => res.iter().map(|x| x.ligand.altloc.clone()).collect::<Vec<String>>(),
        "from_atomn" => res.iter().map(|x| x.ligand.atomn.clone()).collect::<Vec<String>>(),
        "from_atomi" => res.iter().map(|x| x.ligand.atomi as u64).collect::<Vec<u64>>(),
        "to_chain" => res.iter().map(|x| x.receptor.chain.clone()).collect::<Vec<String>>(),
        "to_resn" => res.iter().map(|x| x.receptor.resn.clone()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|x| x.receptor.resi as i64).collect::<Vec<i64>>(),
        "to_insertion" => res.iter().map(|x| x.receptor.insertion.clone()).collect::<Vec<String>>(),
        "to_altloc" => res.iter().map(|x| x.receptor.altloc.clone()).collect::<Vec<String>>(),
        "to_atomn" => res.iter().map(|x| x.receptor.atomn.clone()).collect::<Vec<String>>(),
        "to_atomi" => res.iter().map(|x| x.receptor.atomi as u64).collect::<Vec<u64>>(),
    )
    .unwrap()
}

/// Convert sidechain statistics into a Polars `DataFrame`.
pub(crate) fn sc_results_to_df(
    res: &HashMap<(ResidueId, ResidueId), (f64, f64, f64)>,
) -> DataFrame {
    df!(
        "model" => res.keys().map(|k| k.0.model as u64).collect::<Vec<u64>>(),
        "from_chain" => res.keys().map(|k| k.0.chain.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.keys().map(|k| k.0.resi as i64).collect::<Vec<i64>>(),
        "from_insertion" => res.keys().map(|k| k.0.insertion.to_owned()).collect::<Vec<String>>(),
        "from_altloc" => res.keys().map(|k| k.0.altloc.to_owned()).collect::<Vec<String>>(),
        "to_chain" => res.keys().map(|k| k.1.chain.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.keys().map(|k| k.1.resi as i64).collect::<Vec<i64>>(),
        "to_insertion" => res.keys().map(|k| k.1.insertion.to_owned()).collect::<Vec<String>>(),
        "to_altloc" => res.keys().map(|k| k.1.altloc.to_owned()).collect::<Vec<String>>(),
        "sc_centroid_dist" => res.values().map(|v| v.0 as f32).collect::<Vec<f32>>(),
        "sc_dihedral" => res.values().map(|v| v.1 as f32).collect::<Vec<f32>>(),
        "sc_centroid_angle" => res.values().map(|v| v.2 as f32).collect::<Vec<f32>>(),
    )
    .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    fn default_contacts(pdb: &PDB, groups: &str, vdw_comp: f64, cutoff: f64) -> DataFrame {
        get_contacts(
            pdb,
            groups,
            vdw_comp,
            cutoff,
            ProtonationMode::AllCharged,
            7.4,
        )
        .unwrap()
        .value
    }

    fn two_carbon_atoms(distance: f64) -> PDB {
        let input = format!(
            "ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
             ATOM      2  CB  ALA B   1       {distance:>5.3}   0.000   0.000  1.00 20.00           C  \n\
             END                                                                             \n"
        );
        ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_bytes()))
            .unwrap()
            .0
    }

    #[test]
    fn ringless_proteins_still_produce_atomic_contacts() {
        let input =
            b"ATOM      1  NZ  LYS A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
ATOM      2  OD1 ASP B   1       2.500   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        let contacts = default_contacts(&pdb, "A/B", 0.1, 6.5);

        assert!(contacts.height() > 0);
    }

    #[test]
    fn atomic_contacts_never_cross_model_boundaries() {
        let input = b"MODEL        1\n\
ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CB  ALA B   1      10.000   0.000   0.000  1.00 20.00           C  \n\
ENDMDL\n\
MODEL        2\n\
ATOM      1  CB  ALA A   1      10.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CB  ALA B   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ENDMDL\nEND\n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        assert_eq!(default_contacts(&pdb, "A/B", 0.1, 6.5).height(), 0);
    }

    #[test]
    fn coordinate_neighbors_without_a_peptide_bond_are_compared() {
        let input =
            b"ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CB  ALA A  10       3.600   0.000   0.000  1.00 20.00           C  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        let contacts = default_contacts(&pdb, "/", 0.1, 6.5);
        assert!(
            contacts
                .column("interaction")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .flatten()
                .any(|value| value == "VanDerWaalsContact")
        );
    }

    #[test]
    fn peptide_bonded_neighbors_are_not_compared() {
        let input =
            b"ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  N   ALA A   2       1.330   0.000   0.000  1.00 20.00           N  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        assert_eq!(default_contacts(&pdb, "/", 0.1, 6.5).height(), 0);
    }

    #[test]
    fn valid_ter_prevents_peptide_neighbor_suppression() {
        let input =
            b"ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
TER       2      ALA A   1\n\
ATOM      3  N   ALA A   2       1.330   0.000   0.000  1.00 20.00           N  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();
        let path =
            std::env::temp_dir().join(format!("arpeggia-valid-ter-{}.pdb", std::process::id()));
        std::fs::write(&path, input).unwrap();
        let metadata = crate::read_metadata(&path).unwrap().value;
        std::fs::remove_file(path).unwrap();

        let contacts = get_contacts_with_metadata(
            &pdb,
            &metadata,
            "/",
            0.1,
            6.5,
            ProtonationMode::AllCharged,
            7.4,
        )
        .unwrap()
        .value;
        assert!(contacts.height() > 0);
    }

    #[test]
    fn distinguishes_contact_distance_regions() {
        for (distance, expected) in [
            (1.0, "StericClash"),
            (1.5, "PotentialCovalent"),
            (2.0, "VanDerWaalsClash"),
            (3.6, "VanDerWaalsContact"),
        ] {
            let contacts = default_contacts(&two_carbon_atoms(distance), "A/B", 0.1, 6.5);
            let interactions = contacts.column("interaction").unwrap().str().unwrap();
            assert!(
                interactions.iter().flatten().any(|value| value == expected),
                "expected {expected} at {distance} Å; got {interactions:?}"
            );
        }
    }

    #[test]
    fn explicit_input_connectivity_takes_precedence_over_distance() {
        let input = b"ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CB  ALA B   1       1.500   0.000   0.000  1.00 20.00           C  \n\
CONECT    1    2\nEND                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();
        let path =
            std::env::temp_dir().join(format!("arpeggia-connectivity-{}.pdb", std::process::id()));
        std::fs::write(&path, input).unwrap();
        let metadata = crate::read_metadata(&path).unwrap().value;
        std::fs::remove_file(path).unwrap();

        let contacts = get_contacts_with_metadata(
            &pdb,
            &metadata,
            "A/B",
            0.1,
            6.5,
            ProtonationMode::AllCharged,
            7.4,
        )
        .unwrap()
        .value;
        let interactions = contacts.column("interaction").unwrap().str().unwrap();
        assert!(
            interactions
                .iter()
                .flatten()
                .any(|value| value == "Covalent")
        );
        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "PotentialCovalent")
        );
    }

    #[test]
    fn histidine_protonation_modes_separate_potential_and_explicit_charge() {
        let input =
            b"ATOM      1  CG  HIS A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  OD1 ASP B   1       3.000   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let (mut pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();
        let interactions = |analysis: &crate::Analysis<DataFrame>| {
            analysis
                .value
                .column("interaction")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .flatten()
                .map(str::to_string)
                .collect::<Vec<_>>()
        };

        let all = analyze_contacts(
            &pdb,
            None,
            "A/B",
            0.1,
            6.5,
            ProtonationMode::AllCharged,
            7.4,
        )
        .unwrap();
        assert!(interactions(&all).contains(&"PotentialIonicBond".into()));
        assert!(
            all.warnings
                .iter()
                .any(|warning| warning.code == crate::WarningCode::UnresolvedHistidine)
        );

        let acidic =
            analyze_contacts(&pdb, None, "A/B", 0.1, 6.5, ProtonationMode::Heuristic, 5.0).unwrap();
        assert!(interactions(&acidic).contains(&"PotentialIonicBond".into()));
        let physiological =
            analyze_contacts(&pdb, None, "A/B", 0.1, 6.5, ProtonationMode::Heuristic, 7.4).unwrap();
        assert!(!interactions(&physiological).contains(&"PotentialIonicBond".into()));

        let explicit_only = analyze_contacts(
            &pdb,
            None,
            "A/B",
            0.1,
            6.5,
            ProtonationMode::ExplicitOnly,
            7.4,
        )
        .unwrap();
        assert!(!interactions(&explicit_only).contains(&"IonicBond".into()));
        assert!(!interactions(&explicit_only).contains(&"PotentialIonicBond".into()));

        pdb.atoms_mut()
            .find(|atom| atom.name() == "CG")
            .unwrap()
            .set_charge(1);
        let explicit = analyze_contacts(
            &pdb,
            None,
            "A/B",
            0.1,
            6.5,
            ProtonationMode::ExplicitOnly,
            7.4,
        )
        .unwrap();
        assert!(interactions(&explicit).contains(&"IonicBond".into()));
    }

    #[test]
    fn unrelated_residue_hydrogen_cannot_prove_hydrogen_bond() {
        let input =
            b"ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
ATOM      2  HZ1 LYS A   1       1.000   0.000   0.000  1.00 20.00           H  \n\
ATOM      3  OD1 ASP B   1       2.500   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        let contacts = default_contacts(&pdb, "A/B", 0.1, 6.5);
        let interactions = contacts.column("interaction").unwrap().str().unwrap();

        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "HydrogenBond")
        );
        assert!(
            interactions
                .iter()
                .flatten()
                .any(|value| value == "PolarContact")
        );
        let analysis = analyze_contacts(
            &pdb,
            None,
            "A/B",
            0.1,
            6.5,
            ProtonationMode::AllCharged,
            7.4,
        )
        .unwrap();
        assert!(
            analysis
                .warnings
                .iter()
                .any(|warning| warning.code == crate::WarningCode::MissingDonorHydrogen)
        );
    }

    #[test]
    fn explicit_connectivity_can_associate_a_noncanonical_donor_hydrogen() {
        let input = b"ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
ATOM      2  Q1  LYS A   1       1.000   0.000   0.000  1.00 20.00           H  \n\
ATOM      3  OD1 ASP B   1       2.500   0.000   0.000  1.00 20.00           O  \n\
CONECT    1    2\nEND                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();
        let path = std::env::temp_dir().join(format!(
            "arpeggia-donor-connectivity-{}.pdb",
            std::process::id()
        ));
        std::fs::write(&path, input).unwrap();
        let metadata = crate::read_metadata(&path).unwrap().value;
        std::fs::remove_file(path).unwrap();

        let analysis = analyze_contacts(
            &pdb,
            Some(&metadata),
            "A/B",
            0.1,
            6.5,
            ProtonationMode::AllCharged,
            7.4,
        )
        .unwrap();
        assert!(
            analysis
                .value
                .column("interaction")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .flatten()
                .any(|value| value == "HydrogenBond")
        );
        assert!(
            !analysis
                .warnings
                .iter()
                .any(|warning| warning.code == crate::WarningCode::MissingDonorHydrogen)
        );
    }

    #[test]
    fn unrelated_residue_hydrogen_cannot_prove_weak_hydrogen_bond() {
        let input =
            b"ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  HA  ALA A   1       1.000   0.000   0.000  1.00 20.00           H  \n\
ATOM      3  OD1 ASP B   1       2.500   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        let contacts = default_contacts(&pdb, "A/B", 0.1, 6.5);
        let interactions = contacts.column("interaction").unwrap().str().unwrap();

        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "WeakHydrogenBond")
        );
        assert!(
            interactions
                .iter()
                .flatten()
                .any(|value| value == "WeakPolarContact")
        );
    }
}
