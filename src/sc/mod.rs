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

use crate::structure::{parse_groups, prepare_structure_with_chains};
use pdbtbx::PDB;
use sc_calculator::ScCalculator;
use std::collections::HashSet;

pub use surface_generator::SurfaceCalculatorError;

/// Combined and directional shape-complementarity medians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScResult {
    /// Median across both directional surfaces.
    pub sc: f64,
    /// Median score for group 1 evaluated against group 2.
    pub group1: f64,
    /// Median score for group 2 evaluated against group 1.
    pub group2: f64,
}

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
/// let pdb = load_model(&input_file).unwrap().value;
/// let sc_score = get_sc(&pdb, "H,L/A", 0).unwrap().value;
/// println!("Shape complementarity: {:.3}", sc_score);
/// ```
pub fn get_sc(
    pdb: &PDB,
    groups: &str,
    model_num: usize,
) -> Result<crate::Analysis<f64>, SurfaceCalculatorError> {
    Ok(get_sc_details(pdb, groups, model_num)?.map(|result| result.sc))
}

/// Calculate combined and both directional SC scores.
pub fn get_sc_details(
    pdb: &PDB,
    groups: &str,
    model_num: usize,
) -> Result<crate::Analysis<ScResult>, SurfaceCalculatorError> {
    let models = pdb
        .models()
        .map(|model| model.serial_number())
        .collect::<Vec<_>>();
    if models.is_empty() || (model_num != 0 && !models.contains(&model_num)) {
        return Err(SurfaceCalculatorError::InvalidInput(format!(
            "model {model_num} does not exist"
        )));
    }
    let mut pdb = pdb.clone();
    let warnings = crate::structure::select_conformers(&mut pdb);
    // Get all chains in the PDB
    let all_chains: HashSet<String> = pdb.chains().map(|c| c.id().to_string()).collect();

    // Parse groups
    let (group1_chains, group2_chains) = parse_groups(&all_chains, groups)
        .map_err(|error| SurfaceCalculatorError::InvalidInput(error.to_string()))?;
    if !group1_chains.is_disjoint(&group2_chains) {
        return Err(SurfaceCalculatorError::InvalidInput(
            "SC chain groups must be disjoint".into(),
        ));
    }

    let all_selected_chains: HashSet<String> =
        group1_chains.union(&group2_chains).cloned().collect();
    let pdb_filtered = prepare_structure_with_chains(&pdb, model_num, true, &all_selected_chains);

    // Initialize the calculator with thread settings
    let mut calc = ScCalculator::new();

    // Load atoms from PDB into the calculator
    // Each contains a `molecule_id` that is 0 for group1, 1 for group2
    calc.add_atoms(&pdb_filtered, &group1_chains, &group2_chains)?;

    // Calculate SC
    let results = calc.calc()?;
    Ok(crate::Analysis::new(
        ScResult {
            sc: results.sc,
            group1: results.surfaces[0].s_median,
            group2: results.surfaces[1].s_median,
        },
        warnings,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{load_model, run_with_threads};

    fn load_multi_chain() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/6bft.pdb");
        load_model(&path).unwrap().value
    }

    #[test]
    fn test_h_vs_l() {
        let pdb = load_multi_chain();
        let result = match run_with_threads(0, || get_sc_details(&pdb, "H/L", 0)) {
            Ok(value) => value.value,
            Err(e) => panic!("Error calculating SC: {:?}", e),
        };

        // [WARNING] These corrected-port regression values are not an independent
        // absolute reference: pinned sc-rs reports 0.714361 for this input.
        let expected_sc = 0.725407;
        assert!(
            (result.sc - expected_sc).abs() < 0.000_005,
            "Expected SC around {expected_sc}, got {result:?}",
        );
        assert!((result.group1 - 0.724196).abs() < 0.000_005);
        assert!((result.group2 - 0.726618).abs() < 0.000_005);
    }

    #[test]
    fn sc_rejects_overlapping_groups() {
        let pdb = load_multi_chain();
        assert!(matches!(
            get_sc(&pdb, "H/H", 0),
            Err(SurfaceCalculatorError::InvalidInput(_))
        ));
    }

    #[test]
    fn test_h_vs_c() {
        let pdb = load_multi_chain();
        let sc_value = match run_with_threads(0, || get_sc(&pdb, "H/C", 0)) {
            Ok(value) => value.value,
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
            Ok(value) => value.value,
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
    fn chains_without_interface_return_an_error() {
        let pdb = load_multi_chain();
        let result = run_with_threads(1, || get_sc(&pdb, "H/B", 0));
        assert!(matches!(result, Err(SurfaceCalculatorError::NoInterface)));
    }
}
