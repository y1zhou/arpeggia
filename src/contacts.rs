//! Contact analysis functionality for protein structures.
//!
//! This module provides functions for analyzing atomic and ring contacts
//! in PDB and mmCIF structures.

use crate::interactions::{InteractionComplex, Interactions, ResultEntry};
use crate::residues::ResidueId;
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
/// A Polars DataFrame containing all identified contacts with columns:
/// - model, interaction, distance
/// - from_chain, from_resn, from_resi, from_insertion, from_altloc, from_atomn, from_atomi
/// - to_chain, to_resn, to_resi, to_insertion, to_altloc, to_atomn, to_atomi
/// - sc_centroid_dist, sc_dihedral, sc_centroid_angle
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
    let (i_complex, build_ring_err) =
        InteractionComplex::new(pdb, groups, vdw_comp, dist_cutoff).unwrap();

    if !build_ring_err.is_empty() {
        build_ring_err.iter().for_each(|e| warn!("{e}"));
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
            Default::default(),
        )
        .collect()
        .unwrap()
}

/// Convert a slice of ResultEntry into a Polars DataFrame.
pub(crate) fn results_to_df(res: &[ResultEntry]) -> DataFrame {
    df!(
        "model" => res.iter().map(|x| x.model as u32).collect::<Vec<u32>>(),
        "interaction" => res.iter().map(|x| x.interaction.to_string()).collect::<Vec<String>>(),
        "distance" => res.iter().map(|x| x.distance as f32).collect::<Vec<f32>>(),
        "from_chain" => res.iter().map(|x| x.ligand.chain.to_owned()).collect::<Vec<String>>(),
        "from_resn" => res.iter().map(|x| x.ligand.resn.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.iter().map(|x| x.ligand.resi as i32).collect::<Vec<i32>>(),
        "from_insertion" => res.iter().map(|x| x.ligand.insertion.to_owned()).collect::<Vec<String>>(),
        "from_altloc" => res.iter().map(|x| x.ligand.altloc.to_owned()).collect::<Vec<String>>(),
        "from_atomn" => res.iter().map(|x| x.ligand.atomn.to_owned()).collect::<Vec<String>>(),
        "from_atomi" => res.iter().map(|x| x.ligand.atomi as i32).collect::<Vec<i32>>(),
        "to_chain" => res.iter().map(|x| x.receptor.chain.to_owned()).collect::<Vec<String>>(),
        "to_resn" => res.iter().map(|x| x.receptor.resn.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.iter().map(|x| x.receptor.resi as i32).collect::<Vec<i32>>(),
        "to_insertion" => res.iter().map(|x| x.receptor.insertion.to_owned()).collect::<Vec<String>>(),
        "to_altloc" => res.iter().map(|x| x.receptor.altloc.to_owned()).collect::<Vec<String>>(),
        "to_atomn" => res.iter().map(|x| x.receptor.atomn.to_owned()).collect::<Vec<String>>(),
        "to_atomi" => res.iter().map(|x| x.receptor.atomi as i32).collect::<Vec<i32>>(),
    )
    .unwrap()
}

/// Convert sidechain statistics into a Polars DataFrame.
pub(crate) fn sc_results_to_df(res: &HashMap<(ResidueId, ResidueId), (f64, f64, f64)>) -> DataFrame {
    df!(
        "model" => res.keys().map(|k| k.0.model as u32).collect::<Vec<u32>>(),
        "from_chain" => res.keys().map(|k| k.0.chain.to_owned()).collect::<Vec<String>>(),
        "from_resi" => res.keys().map(|k| k.0.resi as i32).collect::<Vec<i32>>(),
        "from_insertion" => res.keys().map(|k| k.0.insertion.to_owned()).collect::<Vec<String>>(),
        "from_altloc" => res.keys().map(|k| k.0.altloc.to_owned()).collect::<Vec<String>>(),
        "to_chain" => res.keys().map(|k| k.1.chain.to_owned()).collect::<Vec<String>>(),
        "to_resi" => res.keys().map(|k| k.1.resi as i32).collect::<Vec<i32>>(),
        "to_insertion" => res.keys().map(|k| k.1.insertion.to_owned()).collect::<Vec<String>>(),
        "to_altloc" => res.keys().map(|k| k.1.altloc.to_owned()).collect::<Vec<String>>(),
        "sc_centroid_dist" => res.values().map(|v| v.0 as f32).collect::<Vec<f32>>(),
        "sc_dihedral" => res.values().map(|v| v.1 as f32).collect::<Vec<f32>>(),
        "sc_centroid_angle" => res.values().map(|v| v.2 as f32).collect::<Vec<f32>>(),
    )
    .unwrap()
}
