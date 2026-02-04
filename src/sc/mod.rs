//! Shape Complementarity (SC) calculation module.
//!
//! Implements the algorithm from Lawrence & Colman (1993) for calculating
//! shape complementarity between protein interfaces.
//!
//! Based on https://github.com/cytokineking/sc-rs

pub mod atomic_radii;
pub mod sc_calculator;
pub mod settings;
pub mod surface_generator;
pub mod types;
pub mod vector3;

use crate::sasa::{filter_pdb_by_model, prepare_pdb_for_sasa};
use crate::utils::parse_groups;
use pdbtbx::{
    ContainsAtomConformer, ContainsAtomConformerResidue, ContainsAtomConformerResidueChain, PDB,
};
use sc_calculator::ScCalculator;
use std::collections::HashSet;
use types::Atom;
use vector3::Vec3;

pub use surface_generator::SurfaceCalculatorError;

/// Calculate Shape Complementarity (SC) between two chain groups.
///
/// SC measures geometric fit at protein-protein interfaces.
/// Higher SC values indicate better shape complementarity.
/// - Typical antibody-antigen: 0.64-0.68
/// - Typical protein-protein: 0.5-0.7
/// - SC > 0.7 indicates very good complementarity
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `groups` - Chain groups specification (e.g., "H,L/A" or "A/B")
/// * `model_num` - Model number to analyze (0 for first model)
/// * `threads` - Number of threads for parallel calculations (0 for auto)
///
/// # Returns
///
/// The SC score as f64.
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_sc};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
/// let sc_score = get_sc(&pdb, "H,L/A", 0, 0).unwrap();
/// println!("Shape complementarity: {:.3}", sc_score);
/// ```
pub fn get_sc(
    pdb: &PDB,
    groups: &str,
    model_num: usize,
    threads: usize,
) -> Result<f64, SurfaceCalculatorError> {
    // Get all chains in the PDB
    let all_chains: HashSet<String> = pdb.chains().map(|c| c.id().to_string()).collect();

    // Parse groups
    let (group1_chains, group2_chains) = parse_groups(&all_chains, groups);

    // Combine all chains from both groups for filtering
    let all_selected_chains: String = group1_chains
        .iter()
        .chain(group2_chains.iter())
        .cloned()
        .collect::<Vec<_>>()
        .join(",");

    // Prepare PDB: remove hydrogens, solvent, and ions; filter to selected chains
    let pdb_prepared = prepare_pdb_for_sasa(pdb, true, true, &all_selected_chains);

    // Filter to specified model
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    // Initialize the calculator with thread settings
    let mut calc = ScCalculator::new();
    calc.set_threads(threads);

    // Load atoms from PDB into the calculator using filter_map/collect pattern
    // Each entry is (molecule_id, Atom) where molecule_id is 0 for group1, 1 for group2
    let atom_assignments: Vec<Atom> = pdb_filtered
        .atoms_with_hierarchy()
        .filter_map(|hier| {
            let chain_id = hier.chain().id().to_string();

            // Determine which molecule (0 or 1) based on chain group
            let molecule = if group1_chains.contains(&chain_id) {
                0
            } else if group2_chains.contains(&chain_id) {
                1
            } else {
                return None; // Skip atoms not in either group
            };

            let atom = hier.atom();
            let residue = hier.residue();
            let pos = atom.pos();

            Some(Atom {
                natom: 0, // Will be assigned by add_atom
                molecule,
                radius: 0.0, // Will be assigned from radii table
                density: 0.0,
                attention: types::Attention::Buried,
                accessible: false,
                atom: atom.name().to_string(),
                residue: residue.name().unwrap_or("UNK").to_string(),
                coor: Vec3::new(pos.0, pos.1, pos.2),
                neighbor_indices: Vec::new(),
                buried_by_indices: Vec::new(),
                elem_radius: atom.element().unwrap().atomic_radius().van_der_waals,
            })
        })
        .collect();

    // Add atoms to calculator
    for atom in atom_assignments {
        calc.add_atom(atom)?;
    }

    // Calculate SC
    let results = calc.calc()?;

    Ok(results.sc)
}
