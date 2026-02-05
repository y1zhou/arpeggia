//! Shape Complementarity (SC) calculation module.
//!
//! Implements the algorithm from Lawrence & Colman (1993) for calculating
//! shape complementarity between protein interfaces.
//!
//! Based on <https://github.com/cytokineking/sc-rs>

pub mod atomic_radii;
pub mod sc_calculator;
pub mod settings;
pub mod surface_generator;
pub mod types;
pub mod vector3;

use crate::sasa::{filter_pdb_by_model, prepare_pdb_for_sasa};
use crate::utils::parse_groups;
use pdbtbx::PDB;
use sc_calculator::ScCalculator;
use std::collections::HashSet;

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
/// let sc_score = get_sc(&pdb, "H,L/A", 0).unwrap();
/// println!("Shape complementarity: {:.3}", sc_score);
/// ```
pub fn get_sc(pdb: &PDB, groups: &str, model_num: usize) -> Result<f64, SurfaceCalculatorError> {
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

    // Load atoms from PDB into the calculator
    // Each contains a `molecule_id` that is 0 for group1, 1 for group2
    calc.add_atoms(&pdb_filtered, &group1_chains, &group2_chains)?;

    // Calculate SC
    let results = calc.calc()?;
    Ok(results.sc)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::{load_model, run_with_threads};

    fn load_multi_chain() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/6bft.pdb");
        let (pdb, _) = load_model(&path);
        pdb
    }

    #[test]
    fn test_h_vs_l() {
        let pdb = load_multi_chain();
        let sc_value = match run_with_threads(0, || get_sc(&pdb, "H/L", 0)) {
            Ok(value) => value,
            Err(e) => panic!("Error calculating SC: {:?}", e),
        };

        // Check that the result matches our expections
        let expected_sc = 0.714;
        assert!(
            (sc_value - expected_sc).abs() < 0.05,
            "Expected SC around {expected_sc}, got {sc_value}",
        );
    }

    #[test]
    fn test_h_vs_c() {
        let pdb = load_multi_chain();
        let sc_value = match run_with_threads(0, || get_sc(&pdb, "H/C", 0)) {
            Ok(value) => value,
            Err(e) => panic!("Error calculating SC: {:?}", e),
        };

        // Check that the result matches our expections
        let expected_sc = 0.785;
        assert!(
            (sc_value - expected_sc).abs() < 0.05,
            "Expected SC around {expected_sc}, got {sc_value}",
        );
    }

    #[test]
    fn test_hl_vs_cg() {
        let pdb = load_multi_chain();
        let sc_value = match run_with_threads(0, || get_sc(&pdb, "H,L/C,G", 0)) {
            Ok(value) => value,
            Err(e) => panic!("Error calculating SC: {:?}", e),
        };

        // Check that the result matches our expections
        let expected_sc = 0.745;
        assert!(
            (sc_value - expected_sc).abs() < 0.05,
            "Expected SC around {expected_sc}, got {sc_value}",
        );
    }

    #[test]
    #[should_panic(expected = "No molecular dots generated")]
    fn test_chains_without_interface() {
        let pdb = load_multi_chain();
        let _ = match run_with_threads(1, || get_sc(&pdb, "H/B", 0)) {
            Ok(value) => value,
            Err(e) => panic!("Error calculating SC: {:?}", e),
        };
    }
}
