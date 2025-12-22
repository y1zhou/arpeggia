#![warn(missing_docs)]
#![doc = include_str!("../README.md")]

//! # Arpeggia Library
//!
//! This library provides functionality for analyzing protein-protein interactions
//! in PDB and mmCIF files. It can identify various types of interactions including
//! hydrogen bonds, ionic interactions, aromatic interactions, and more.
//!
//! The library returns results as Polars DataFrames, which can be easily converted
//! to various output formats or used directly in Python via PyO3 bindings.

mod chains;
mod interactions;
mod residues;
mod utils;

// Re-export key public types
pub use interactions::{InteractingEntity, Interaction, Interactions, ResultEntry};
pub use residues::{Plane, ResidueExt, ResidueId};
pub use utils::{load_model, parse_groups, write_df_to_file, DataFrameFileType};

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
        interactions::InteractionComplex::new(pdb, groups, vdw_comp, dist_cutoff).unwrap();

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

/// Calculate solvent accessible surface area (SASA) for each atom in a PDB structure.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
///
/// # Returns
///
/// A Polars DataFrame with columns:
/// - atomi, sasa
/// - chain, resn, resi, insertion, altloc, atomn
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_atom_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
/// let sasa_df = get_atom_sasa(&pdb, 1.4, 100, 0);
/// println!("Calculated SASA for {} atoms", sasa_df.height());
/// ```
pub fn get_atom_sasa(pdb: &PDB, probe_radius: f32, n_points: usize, model_num: usize) -> DataFrame {
    use crate::residues::ResidueExt;
    use rust_sasa::{calculate_sasa_internal, Atom as SASAAtom};

    // If model_num is 0, we use the first model; otherwise use the specified model
    let model_num = if model_num == 0 {
        pdb.models()
            .collect::<Vec<_>>()
            .first()
            .map_or(0, |m| m.serial_number())
    } else {
        model_num
    };

    // Calculate the SASA for each atom
    let atoms = pdb
        .atoms_with_hierarchy()
        .filter(|x| {
            let resn = x.residue().resn().unwrap();
            resn != "O" && resn != "X" && x.model().serial_number() == model_num
        })
        .map(|x| SASAAtom {
            position: nalgebra::Point3::new(
                x.atom().pos().0 as f32,
                x.atom().pos().1 as f32,
                x.atom().pos().2 as f32,
            ),
            radius: x
                .atom()
                .element()
                .unwrap()
                .atomic_radius()
                .van_der_waals
                .unwrap() as f32,
            id: x.atom().serial_number(),
            parent_id: None,
        })
        .collect::<Vec<_>>();
    let atom_sasa = calculate_sasa_internal(&atoms, probe_radius, n_points, true);

    // Create a DataFrame with the results
    let atom_annotations = pdb
        .atoms_with_hierarchy()
        .map(|x| InteractingEntity::from_hier(&x))
        .collect::<Vec<InteractingEntity>>();
    let atom_annot_df = df!(
        "chain" => atom_annotations.iter().map(|x| x.chain.to_owned()).collect::<Vec<String>>(),
        "resn" => atom_annotations.iter().map(|x| x.resn.to_owned()).collect::<Vec<String>>(),
        "resi" => atom_annotations.iter().map(|x| x.resi as i32).collect::<Vec<i32>>(),
        "insertion" => atom_annotations.iter().map(|x| x.insertion.to_owned()).collect::<Vec<String>>(),
        "altloc" => atom_annotations.iter().map(|x| x.altloc.to_owned()).collect::<Vec<String>>(),
        "atomn" => atom_annotations.iter().map(|x| x.atomn.to_owned()).collect::<Vec<String>>(),
        "atomi" => atom_annotations.iter().map(|x| x.atomi as i32).collect::<Vec<i32>>(),
    )
    .unwrap();

    df!(
        "atomi" => atoms.iter().map(|x| x.id as i32).collect::<Vec<i32>>(),
        "sasa" => atom_sasa
    )
    .unwrap()
    .join(
        &atom_annot_df,
        ["atomi"],
        ["atomi"],
        JoinArgs::new(JoinType::Inner),
        None,
    )
    .unwrap()
    .sort(["chain", "resi", "altloc", "atomi"], Default::default())
    .unwrap()
}

/// Get sequences of all chains in a PDB structure.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
///
/// # Returns
///
/// A HashMap mapping chain IDs to their sequences as strings.
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_sequences};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
/// let sequences = get_sequences(&pdb);
/// for (chain_id, seq) in sequences {
///     println!("Chain {}: {}", chain_id, seq);
/// }
/// ```
pub fn get_sequences(pdb: &PDB) -> HashMap<String, String> {
    use crate::chains::ChainExt;

    pdb.chains()
        .map(|chain| (chain.id().to_string(), chain.pdb_seq().join("")))
        .collect()
}

// Helper functions (kept private)

fn results_to_df(res: &[ResultEntry]) -> DataFrame {
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

fn sc_results_to_df(res: &HashMap<(ResidueId, ResidueId), (f64, f64, f64)>) -> DataFrame {
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

// Python bindings module (only compiled when python feature is enabled)
#[cfg(feature = "python")]
mod python;

#[cfg(feature = "python")]
pub use python::*;
