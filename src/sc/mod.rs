//! Shape Complementarity (SC) calculation module.
//!
//! Implements the algorithm from Lawrence & Colman (1993) for calculating
//! shape complementarity between protein interfaces.
//!
//! Based on <https://github.com/cytokineking/sc-rs>

pub mod atomic_radii;
pub mod sc_calculator;
pub mod surface_generator;
pub mod types;
pub mod vector3;

use crate::structure::{parse_groups, prepare_structure_with_chains};
use pdbtbx::PDB;
use sc_calculator::ScCalculator;
use std::collections::HashSet;

const PROBE_RADIUS: f64 = 1.7;
const DOT_DENSITY: f64 = 15.0;
const PERIPHERAL_BAND: f64 = 1.5;
const SEPARATION_CUTOFF: f64 = 8.0;
const GAUSSIAN_WEIGHT: f64 = 0.5;

pub use surface_generator::SurfaceCalculatorError;

impl From<SurfaceCalculatorError> for crate::ArpeggiaError {
    fn from(error: SurfaceCalculatorError) -> Self {
        match error {
            SurfaceCalculatorError::InvalidInput(message) => Self::InvalidArgument(message),
            error => Self::Calculation(error.to_string()),
        }
    }
}

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
) -> crate::ArpeggiaResult<crate::Analysis<f64>> {
    Ok(get_sc_details(pdb, groups, model_num)?.map(|result| result.sc))
}

/// Calculate combined and both directional SC scores.
pub fn get_sc_details(
    pdb: &PDB,
    groups: &str,
    model_num: usize,
) -> crate::ArpeggiaResult<crate::Analysis<ScResult>> {
    let selected_model = crate::structure::selected_model(pdb, model_num)?;
    let mut pdb = pdb.clone();
    let warnings = crate::structure::select_conformers(&mut pdb);
    // Get all chains in the PDB
    let all_chains: HashSet<String> = selected_model
        .chains()
        .map(|chain| chain.id().to_string())
        .collect();

    // Parse groups
    let (group1_chains, group2_chains) = parse_groups(&all_chains, groups)?;
    if !group1_chains.is_disjoint(&group2_chains) {
        return Err(crate::ArpeggiaError::InvalidArgument(
            "SC chain groups must be disjoint".into(),
        ));
    }

    let all_selected_chains: HashSet<String> =
        group1_chains.union(&group2_chains).cloned().collect();
    let pdb_filtered = prepare_structure_with_chains(&pdb, model_num, true, &all_selected_chains);

    let mut calc = ScCalculator::default();

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

        // The pinned sc-rs core reproduces these values when given the same atom
        // population. Its stock CLI reports 0.714361 because its name-based
        // hydrogen filter also removes the oxygen atom named `OH` in tyrosine.
        let expected_sc = 0.725407;
        assert!(
            (result.sc - expected_sc).abs() < 0.000_005,
            "Expected SC around {expected_sc}, got {result:?}",
        );
        assert!((result.group1 - 0.724196).abs() < 0.000_005);
        assert!((result.group2 - 0.726618).abs() < 0.000_005);
    }

    #[test]
    fn histidine_aliases_use_canonical_sc_radii() {
        let root = env!("CARGO_MANIFEST_DIR");
        let input = std::fs::read_to_string(format!("{root}/test-data/6bft.pdb")).unwrap();
        let alias_input = input.replace(" HIS L", " HID L");
        let canonical = load_model(&format!("{root}/test-data/6bft.pdb"))
            .unwrap()
            .value;
        let alias = pdbtbx::ReadOptions::default()
            .set_format(pdbtbx::Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(alias_input.as_bytes()))
            .unwrap()
            .0;

        let canonical = get_sc_details(&canonical, "H/L", 0).unwrap().value;
        let alias = get_sc_details(&alias, "H/L", 0).unwrap().value;
        assert_eq!(alias, canonical);
    }

    #[test]
    fn sc_rejects_overlapping_groups() {
        let pdb = load_multi_chain();
        assert!(matches!(
            get_sc(&pdb, "H/H", 0),
            Err(crate::ArpeggiaError::InvalidArgument(_))
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
        assert!(matches!(result, Err(crate::ArpeggiaError::Calculation(_))));
    }

    #[test]
    fn trimmed_empty_interface_returns_an_error() {
        let input =
            b"ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CB  ALA B   1       7.500   0.000   0.000  1.00 20.00           C  \n\
END                                                                             \n";
        let pdb = pdbtbx::ReadOptions::default()
            .set_format(pdbtbx::Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(input.as_slice()))
            .unwrap()
            .0;

        assert!(matches!(
            get_sc(&pdb, "A/B", 0),
            Err(crate::ArpeggiaError::Calculation(message))
                if message.contains("interface")
        ));
    }

    #[test]
    fn missing_radius_is_a_public_calculation_error() {
        let input =
            b"ATOM      1  QQ  ALA A   1       0.000   0.000   0.000  1.00 20.00          RN  \n\
ATOM      2  CB  ALA B   1       3.000   0.000   0.000  1.00 20.00           C  \n\
END                                                                             \n";
        let pdb = pdbtbx::ReadOptions::default()
            .set_format(pdbtbx::Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(input.as_slice()))
            .unwrap()
            .0;

        assert!(matches!(
            get_sc(&pdb, "A/B", 0),
            Err(crate::ArpeggiaError::Calculation(message))
                if message.contains("van der Waals radius")
        ));
    }

    #[test]
    fn sc_validates_chains_in_the_selected_model() {
        let input = b"MODEL        7\n\
ATOM      1  CB  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CB  ALA C   1       3.000   0.000   0.000  1.00 20.00           C  \n\
ENDMDL\n\
MODEL        9\n\
ATOM      1  CB  ALA B   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CB  ALA D   1       3.000   0.000   0.000  1.00 20.00           C  \n\
ENDMDL\nEND\n";
        let pdb = pdbtbx::ReadOptions::default()
            .set_format(pdbtbx::Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(std::io::BufReader::new(input.as_slice()))
            .unwrap()
            .0;

        assert!(matches!(
            get_sc(&pdb, "B/D", 7),
            Err(crate::ArpeggiaError::InvalidArgument(_))
        ));
    }
}
