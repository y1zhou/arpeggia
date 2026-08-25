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
pub use hydrophobic::find_hydrophobic_contact;
pub use ionic::{find_ionic_bond_with_protonation, find_ionic_repulsion_with_protonation};
pub use structs::*;
pub use vdw::find_vdw_contact;

use pdbtbx::*;
use polars::prelude::*;
use tracing::debug;

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
    warnings.extend(protonation_warnings(
        &prepared,
        protonation,
        &complex.ligand,
        &complex.receptor,
    ));
    let missing_donors = hbond::count_donors_without_explicit_hydrogen(
        &prepared,
        complex.resolved_bonds(),
        &complex.ligand,
        &complex.receptor,
    );
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
        contacts_from_complex(complex)?,
        warnings,
    ))
}

fn protonation_warnings(
    pdb: &PDB,
    mode: ProtonationMode,
    ligand: &std::collections::HashSet<String>,
    receptor: &std::collections::HashSet<String>,
) -> Vec<crate::AnalysisWarning> {
    pdb.chains()
        .filter(|chain| ligand.contains(chain.id()) || receptor.contains(chain.id()))
        .flat_map(|chain| chain.residues())
        .filter(|residue| ionic::is_histidine(residue.name().unwrap_or("")))
        .filter_map(|residue| {
            let issue = ionic::histidine_preparation_issue(residue)?;
            let (code, message) = match issue {
                ionic::HistidinePreparationIssue::Unresolved => (
                    crate::WarningCode::UnresolvedHistidine,
                    format!(
                        "histidine {} has no explicit protonation evidence; policy {mode:?} applied",
                        residue.serial_number()
                    ),
                ),
                ionic::HistidinePreparationIssue::Inconsistent => (
                    crate::WarningCode::InconsistentHistidine,
                    format!(
                        "histidine {} residue name and ring hydrogens disagree",
                        residue.serial_number()
                    ),
                ),
            };
            Some(crate::AnalysisWarning::new(code, message))
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

fn contacts_from_complex(i_complex: InteractionComplex<'_>) -> crate::ArpeggiaResult<DataFrame> {
    // Information on the sequence of the chains in the model
    debug!(
        "Parsed ligand chains {lig:?}; receptor chains {receptor:?}",
        lig = i_complex.ligand,
        receptor = i_complex.receptor
    );

    // Find interactions
    let atomic_contacts = i_complex.get_atomic_contacts();
    let df_atomic = results_to_df(&atomic_contacts)?;
    debug!(
        "Found {} atom-atom contacts\n{}",
        df_atomic.height(),
        df_atomic
    );

    let mut ring_contacts: Vec<ResultEntry> = Vec::new();
    ring_contacts.extend(i_complex.get_ring_atom_contacts());
    ring_contacts.extend(i_complex.get_ring_ring_contacts());
    let df_ring = results_to_df(&ring_contacts)?;
    debug!("Found {} ring contacts\n{}", df_ring.height(), df_ring);

    // Annotate sidechain centroid distances and dihedrals
    let sc_dist_dihedrals = i_complex.collect_sc_stats(&atomic_contacts, &ring_contacts);

    let df_sc_stats = sc_results_to_df(&sc_dist_dihedrals)?;

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
    Ok(df_atomic
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
        .unwrap())
}

/// Convert a slice of `ResultEntry` into a Polars `DataFrame`.
pub(crate) fn results_to_df(res: &[ResultEntry]) -> crate::ArpeggiaResult<DataFrame> {
    Ok(df!(
        "model" => res.iter().map(|x| compact_model_id(x.model)).collect::<crate::ArpeggiaResult<Vec<_>>>()?,
        "interaction" => res.iter().map(|x| x.interaction.to_string()).collect::<Vec<String>>(),
        "distance" => res.iter().map(|x| x.distance as f32).collect::<Vec<f32>>(),
        "from_chain" => res.iter().map(|x| x.ligand.chain.clone()).collect::<Vec<String>>(),
        "from_resn" => res.iter().map(|x| x.ligand.resn.clone()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|x| compact_residue_id(x.ligand.resi)).collect::<crate::ArpeggiaResult<Vec<_>>>()?,
        "from_insertion" => res.iter().map(|x| x.ligand.insertion.clone()).collect::<Vec<String>>(),
        "from_altloc" => res.iter().map(|x| x.ligand.altloc.clone()).collect::<Vec<String>>(),
        "from_atomn" => res.iter().map(|x| x.ligand.atomn.clone()).collect::<Vec<String>>(),
        "from_atomi" => res.iter().map(|x| compact_atom_id(x.ligand.atomi)).collect::<crate::ArpeggiaResult<Vec<_>>>()?,
        "to_chain" => res.iter().map(|x| x.receptor.chain.clone()).collect::<Vec<String>>(),
        "to_resn" => res.iter().map(|x| x.receptor.resn.clone()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|x| compact_residue_id(x.receptor.resi)).collect::<crate::ArpeggiaResult<Vec<_>>>()?,
        "to_insertion" => res.iter().map(|x| x.receptor.insertion.clone()).collect::<Vec<String>>(),
        "to_altloc" => res.iter().map(|x| x.receptor.altloc.clone()).collect::<Vec<String>>(),
        "to_atomn" => res.iter().map(|x| x.receptor.atomn.clone()).collect::<Vec<String>>(),
        "to_atomi" => res.iter().map(|x| compact_atom_id(x.receptor.atomi)).collect::<crate::ArpeggiaResult<Vec<_>>>()?,
    )
    .unwrap())
}

/// Convert sidechain statistics into a Polars `DataFrame`.
pub(crate) fn sc_results_to_df(res: &[SidechainStat<'_>]) -> crate::ArpeggiaResult<DataFrame> {
    Ok(df!(
        "model" => res.iter().map(|(k, _)| compact_model_id(k.0.model)).collect::<crate::ArpeggiaResult<Vec<_>>>()?,
        "from_chain" => res.iter().map(|(k, _)| k.0.chain.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|(k, _)| compact_residue_id(k.0.resi)).collect::<crate::ArpeggiaResult<Vec<_>>>()?,
        "from_insertion" => res.iter().map(|(k, _)| k.0.insertion.to_owned()).collect::<Vec<String>>(),
        "from_altloc" => res.iter().map(|(k, _)| k.0.altloc.to_owned()).collect::<Vec<String>>(),
        "to_chain" => res.iter().map(|(k, _)| k.1.chain.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|(k, _)| compact_residue_id(k.1.resi)).collect::<crate::ArpeggiaResult<Vec<_>>>()?,
        "to_insertion" => res.iter().map(|(k, _)| k.1.insertion.to_owned()).collect::<Vec<String>>(),
        "to_altloc" => res.iter().map(|(k, _)| k.1.altloc.to_owned()).collect::<Vec<String>>(),
        "sc_centroid_dist" => res.iter().map(|(_, v)| v.0 as f32).collect::<Vec<f32>>(),
        "sc_dihedral" => res.iter().map(|(_, v)| v.1 as f32).collect::<Vec<f32>>(),
        "sc_centroid_angle" => res.iter().map(|(_, v)| v.2 as f32).collect::<Vec<f32>>(),
    )
    .unwrap())
}

fn compact_model_id(value: usize) -> crate::ArpeggiaResult<u32> {
    u32::try_from(value)
        .map_err(|_| crate::ArpeggiaError::Calculation("model identifier exceeds UInt32".into()))
}

fn compact_residue_id(value: isize) -> crate::ArpeggiaResult<i32> {
    i32::try_from(value)
        .map_err(|_| crate::ArpeggiaError::Calculation("residue identifier exceeds Int32".into()))
}

fn compact_atom_id(value: usize) -> crate::ArpeggiaResult<u32> {
    u32::try_from(value)
        .map_err(|_| crate::ArpeggiaError::Calculation("atom identifier exceeds UInt32".into()))
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
    fn contact_identifier_columns_use_compact_dtypes() {
        let contacts = default_contacts(&two_carbon_atoms(3.6), "A/B", 0.1, 6.5);

        assert_eq!(contacts.column("model").unwrap().dtype(), &DataType::UInt32);
        for column in ["from_resi", "to_resi"] {
            assert_eq!(contacts.column(column).unwrap().dtype(), &DataType::Int32);
        }
        for column in ["from_atomi", "to_atomi"] {
            assert_eq!(contacts.column(column).unwrap().dtype(), &DataType::UInt32);
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
        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "PotentialDisulfide")
        );
        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "Disulfide")
        );
    }

    #[test]
    fn prepared_2eab_ssbond_precedes_conect() {
        let input = include_str!("../../test-data/2eab-disulfide.pdb")
            .replace("END\n", "CONECT 8893 9228\nEND\n");
        let path = std::env::temp_dir().join(format!(
            "arpeggia-ssbond-precedence-{}.pdb",
            std::process::id()
        ));
        std::fs::write(&path, input).unwrap();
        let pdb = crate::load_model(path.to_str().unwrap()).unwrap().value;
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
        let interactions = contacts.column("interaction").unwrap().str().unwrap();
        assert!(
            interactions
                .iter()
                .flatten()
                .any(|value| value == "Disulfide")
        );
        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "Covalent")
        );
        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "PotentialDisulfide")
        );
        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "PotentialCovalent")
        );
    }

    #[test]
    fn undeclared_disulfide_geometry_is_potential_disulfide() {
        let input = include_str!("../../test-data/2eab-disulfide.pdb")
            .lines()
            .filter(|line| !line.starts_with("SSBOND"))
            .collect::<Vec<_>>()
            .join("\n");
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_bytes()))
            .unwrap();

        let contacts = default_contacts(&pdb, "/", 0.1, 6.5);
        let interactions = contacts.column("interaction").unwrap().str().unwrap();
        assert!(
            interactions
                .iter()
                .flatten()
                .any(|value| value == "PotentialDisulfide"),
            "expected the undeclared CYS SG pair to be PotentialDisulfide; got {interactions:?}"
        );
    }

    #[test]
    fn undeclared_cysteines_outside_disulfide_geometry_remain_potential_covalent() {
        let input =
            b"ATOM      1  CB  CYS A   1       0.000   1.000   0.000  1.00 20.00           C  \n\
ATOM      2  SG  CYS A   1       0.000   0.000   0.000  1.00 20.00           S  \n\
ATOM      3  SG  CYS B   1       2.030   0.000   0.000  1.00 20.00           S  \n\
ATOM      4  CB  CYS B   1       2.030   1.000   0.000  1.00 20.00           C  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        let contacts = default_contacts(&pdb, "A/B", 0.1, 6.5);
        let interactions = contacts.column("interaction").unwrap().str().unwrap();
        assert!(
            interactions
                .iter()
                .flatten()
                .any(|value| value == "PotentialCovalent")
        );
        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "PotentialDisulfide")
        );
    }

    #[test]
    fn discarded_altloc_connectivity_does_not_apply_to_selected_atoms() {
        let input =
            b"ATOM      1  SG ACYS A   1       0.000   0.000   0.000  0.60 20.00           S  \n\
ATOM      2  SG BCYS A   1      20.000   0.000   0.000  0.40 20.00           S  \n\
ATOM      3  SG ACYS B   2       2.030   0.000   0.000  0.60 20.00           S  \n\
ATOM      4  SG BCYS B   2      22.030   0.000   0.000  0.40 20.00           S  \n\
LINK         SG BCYS A   1                 SG BCYS B   2     1555   1555  2.03\n\
END                                                                             \n";
        let mut pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap()
            .0;
        crate::structure::select_conformers(&mut pdb);
        let path = std::env::temp_dir().join(format!(
            "arpeggia-altloc-connectivity-{}.pdb",
            std::process::id()
        ));
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
                .any(|value| value == "PotentialCovalent")
        );
        assert!(
            !interactions
                .iter()
                .flatten()
                .any(|value| value == "Covalent")
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
    fn all_charged_overrides_neutral_histidine_aliases() {
        let input =
            b"ATOM      1  ND1 HID A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
ATOM      2  OD1 ASP B   1       3.800   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap()
            .0;
        let interactions = default_contacts(&pdb, "A/B", 0.1, 6.5);

        assert!(
            interactions
                .column("interaction")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .flatten()
                .any(|value| value == "PotentialIonicBond")
        );
    }

    #[test]
    fn inconsistent_histidine_alias_emits_a_warning() {
        let input =
            b"ATOM      1  ND1 HID A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
ATOM      2  HE2 HID A   1       1.000   0.000   0.000  1.00 20.00           H  \n\
ATOM      3  OD1 ASP B   1       3.800   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap()
            .0;
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
                .any(|warning| warning.code == crate::WarningCode::InconsistentHistidine)
        );
    }

    #[test]
    fn excluded_chains_do_not_emit_contact_preparation_warnings() {
        let input =
            b"ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CB  VAL B   1       3.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      3  ND1 HIS C   1       6.000   0.000   0.000  1.00 20.00           N  \n\
END                                                                             \n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap()
            .0;
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

        assert!(!analysis.warnings.iter().any(|warning| {
            matches!(
                warning.code,
                crate::WarningCode::UnresolvedHistidine | crate::WarningCode::MissingDonorHydrogen
            )
        }));
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
    fn histidine_tautomers_have_specific_donors_and_acceptors() {
        for (residue, donor, hydrogen) in [
            ("HID", "ND1", "HD1"),
            ("HSD", "ND1", "HD1"),
            ("HIE", "NE2", "HE2"),
            ("HSE", "NE2", "HE2"),
            ("HIP", "ND1", "HD1"),
            ("HSP", "ND1", "HD1"),
        ] {
            let donor_input = format!(
                "ATOM      1 {donor:>4} {residue} A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
                 ATOM      2 {hydrogen:>4} {residue} A   1       1.000   0.000   0.000  1.00 20.00           H  \n\
                 ATOM      3  OD1 ASN B   1       2.500   0.000   0.000  1.00 20.00           O  \n\
                 END                                                                             \n"
            );
            let (pdb, _) = ReadOptions::default()
                .set_format(Format::Pdb)
                .set_only_atomic_coords(true)
                .read_raw(BufReader::new(donor_input.as_bytes()))
                .unwrap();
            let interactions = default_contacts(&pdb, "A/B", 0.1, 6.5);
            assert!(
                interactions
                    .column("interaction")
                    .unwrap()
                    .str()
                    .unwrap()
                    .iter()
                    .flatten()
                    .any(|value| value == "HydrogenBond"),
                "{residue} {donor} should donate"
            );

            let acceptor_input = format!(
                "ATOM      1  NE2 GLN A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
                 ATOM      2 HE21 GLN A   1       1.000   0.000   0.000  1.00 20.00           H  \n\
                 ATOM      3 {donor:>4} {residue} B   1       2.500   0.000   0.000  1.00 20.00           N  \n\
                 END                                                                             \n"
            );
            let (pdb, _) = ReadOptions::default()
                .set_format(Format::Pdb)
                .set_only_atomic_coords(true)
                .read_raw(BufReader::new(acceptor_input.as_bytes()))
                .unwrap();
            let interactions = default_contacts(&pdb, "A/B", 0.1, 6.5);
            assert!(
                !interactions
                    .column("interaction")
                    .unwrap()
                    .str()
                    .unwrap()
                    .iter()
                    .flatten()
                    .any(|value| value == "HydrogenBond"),
                "{residue} {donor} should not accept"
            );
        }
    }

    #[test]
    fn explicit_his_ring_hydrogen_selects_the_tautomer() {
        let input =
            b"ATOM      1  NZ  LYS A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
ATOM      2  HZ1 LYS A   1       1.000   0.000   0.000  1.00 20.00           H  \n\
ATOM      3  ND1 HIS B   1       2.500   0.000   0.000  1.00 20.00           N  \n\
ATOM      4  HD1 HIS B   1       3.500   0.000   0.000  1.00 20.00           H  \n\
END                                                                             \n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();
        let interactions = default_contacts(&pdb, "A/B", 0.1, 6.5);
        assert!(
            !interactions
                .column("interaction")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .flatten()
                .any(|value| value == "HydrogenBond")
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

    #[test]
    fn carbon_without_hydrogen_is_not_a_weak_donor() {
        let input =
            b"ATOM      1  CG  ASN A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  OD1 ASP B   1       3.000   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap()
            .0;
        let interactions = default_contacts(&pdb, "A/B", 0.1, 6.5);

        assert!(
            !interactions
                .column("interaction")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .flatten()
                .any(|value| matches!(value, "WeakHydrogenBond" | "WeakPolarContact"))
        );
    }

    #[test]
    fn terminal_proline_with_explicit_hydrogen_can_donate() {
        let input =
            b"ATOM      1  N   PRO A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
ATOM      2  H1  PRO A   1       1.000   0.000   0.000  1.00 20.00           H  \n\
ATOM      3  OD1 ASP B   1       2.500   0.000   0.000  1.00 20.00           O  \n\
END                                                                             \n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap()
            .0;
        let interactions = default_contacts(&pdb, "A/B", 0.1, 6.5);

        assert!(
            interactions
                .column("interaction")
                .unwrap()
                .str()
                .unwrap()
                .iter()
                .flatten()
                .any(|value| value == "HydrogenBond")
        );
    }
}
