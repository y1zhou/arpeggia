//! Contact analysis functionality for protein structures.
//!
//! This module provides functions for analyzing atomic and ring contacts
//! in PDB and mmCIF structures.
pub mod aromatic;
pub mod chains;
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
pub use ionic::{find_ionic_bond, find_ionic_repulsion};
pub use residues::ResidueId;
pub use structs::*;
pub use vdw::find_vdw_contact;

use pdbtbx::*;
use polars::prelude::*;
use std::collections::HashMap;
use tracing::{debug, warn};

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
/// # Panics
///
/// This function will panic if the `InteractionComplex` cannot be created from the provided PDB structure.
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_contacts};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
/// let contacts_df = get_contacts(&pdb, "/", 0.1, 6.5);
/// println!("Found {} contacts", contacts_df.height());
/// ```
pub fn get_contacts(pdb: &PDB, groups: &str, vdw_comp: f64, dist_cutoff: f64) -> DataFrame {
    let (i_complex, build_ring_err) = InteractionComplex::new(pdb, groups, vdw_comp, dist_cutoff);

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
        col("model"),
        col("from_chain"),
        col("from_resi"),
        col("from_insertion"),
        col("from_altloc"),
        col("to_chain"),
        col("to_resi"),
        col("to_insertion"),
        col("to_altloc"),
    ];
    concat([df_atomic.lazy(), df_ring.lazy()], UnionArgs::default())
        .unwrap()
        // Annotate sidechain stats
        .join(
            df_sc_stats.lazy(),
            &contacts_sc_join_cols,
            &contacts_sc_join_cols,
            JoinArgs::new(JoinType::Left),
        )
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
        .collect()
        .unwrap()
}

/// Convert a slice of `ResultEntry` into a Polars `DataFrame`.
pub(crate) fn results_to_df(res: &[ResultEntry]) -> DataFrame {
    df!(
        "model" => res.iter().map(|x| {
            if let Ok(v) = u32::try_from(x.model) {v} else {
                panic!("Model number {} exceeds u32 limits", x.model);
            }
        }).collect::<Vec<u32>>(),
        "interaction" => res.iter().map(|x| x.interaction.to_string()).collect::<Vec<String>>(),
        "distance" => res.iter().map(|x| x.distance as f32).collect::<Vec<f32>>(),
        "from_chain" => res.iter().map(|x| x.ligand.chain.clone()).collect::<Vec<String>>(),
        "from_resn" => res.iter().map(|x| x.ligand.resn.clone()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|x| {
            if let Ok(v) = i32::try_from(x.ligand.resi) {v} else {
                panic!("Residue index {} exceeds i32 limits", x.ligand.resi);
            }
        }).collect::<Vec<i32>>(),
        "from_insertion" => res.iter().map(|x| x.ligand.insertion.clone()).collect::<Vec<String>>(),
        "from_altloc" => res.iter().map(|x| x.ligand.altloc.clone()).collect::<Vec<String>>(),
        "from_atomn" => res.iter().map(|x| x.ligand.atomn.clone()).collect::<Vec<String>>(),
        "from_atomi" => res.iter().map(|x| {
            if let Ok(v) = i32::try_from(x.ligand.atomi) {v} else {
                panic!("Atom index {} exceeds i32 limits", x.ligand.atomi);
            }
        }).collect::<Vec<i32>>(),
        "to_chain" => res.iter().map(|x| x.receptor.chain.clone()).collect::<Vec<String>>(),
        "to_resn" => res.iter().map(|x| x.receptor.resn.clone()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|x| {
            if let Ok(v) = i32::try_from(x.receptor.resi) {v} else {
                panic!("Residue index {} exceeds i32 limits", x.receptor.resi);
            }
        }).collect::<Vec<i32>>(),
        "to_insertion" => res.iter().map(|x| x.receptor.insertion.clone()).collect::<Vec<String>>(),
        "to_altloc" => res.iter().map(|x| x.receptor.altloc.clone()).collect::<Vec<String>>(),
        "to_atomn" => res.iter().map(|x| x.receptor.atomn.clone()).collect::<Vec<String>>(),
        "to_atomi" => res.iter().map(|x| {
            if let Ok(v) = i32::try_from(x.receptor.atomi) {v} else {
                panic!("Atom index {} exceeds i32 limits", x.receptor.atomi);
            }
        }).collect::<Vec<i32>>(),
    )
    .unwrap()
}

/// Convert sidechain statistics into a Polars `DataFrame`.
pub(crate) fn sc_results_to_df(
    res: &HashMap<(ResidueId, ResidueId), (f64, f64, f64)>,
) -> DataFrame {
    df!(
        "model" => res.keys().map(|k| {
            if let Ok(v) = u32::try_from(k.0.model) {v} else {
                panic!("Model number {} exceeds u32 limits", k.0.model);
            }
        }).collect::<Vec<u32>>(),
        "from_chain" => res.keys().map(|k| k.0.chain.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.keys().map(|k| {
            if let Ok(v) = i32::try_from(k.0.resi) {v} else {
                panic!("Residue index {} exceeds i32 limits", k.0.resi);
            }
        }).collect::<Vec<i32>>(),
        "from_insertion" => res.keys().map(|k| k.0.insertion.to_owned()).collect::<Vec<String>>(),
        "from_altloc" => res.keys().map(|k| k.0.altloc.to_owned()).collect::<Vec<String>>(),
        "to_chain" => res.keys().map(|k| k.1.chain.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.keys().map(|k| {
            if let Ok(v) = i32::try_from(k.1.resi) {v} else {
                panic!("Residue index {} exceeds i32 limits", k.1.resi);
            }
        }).collect::<Vec<i32>>(),
        "to_insertion" => res.keys().map(|k| k.1.insertion.to_owned()).collect::<Vec<String>>(),
        "to_altloc" => res.keys().map(|k| k.1.altloc.to_owned()).collect::<Vec<String>>(),
        "sc_centroid_dist" => res.values().map(|v| v.0 as f32).collect::<Vec<f32>>(),
        "sc_dihedral" => res.values().map(|v| v.1 as f32).collect::<Vec<f32>>(),
        "sc_centroid_angle" => res.values().map(|v| v.2 as f32).collect::<Vec<f32>>(),
    )
    .unwrap()
}
