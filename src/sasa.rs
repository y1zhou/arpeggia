//! Solvent Accessible Surface Area (SASA) calculations.
//!
//! This module provides functions for calculating SASA at different levels
//! of granularity (atom, residue, chain) and related metrics like dSASA
//! (buried surface area) and relative SASA.

use crate::interactions::InteractingEntity;
use crate::utils::{parse_groups, sum_sasa};
use pdbtbx::*;
use polars::prelude::*;
use std::collections::HashSet;

/// Filter a PDB structure to keep only the specified model.
///
/// Creates a clone of the PDB and removes all models except the one at the specified index.
/// If model_num is 0, the first model is used. Otherwise, the model with the matching
/// serial number is kept.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `model_num` - Model number to keep (0 for first model)
///
/// # Returns
///
/// A filtered PDB structure containing only the specified model.
pub(crate) fn filter_pdb_by_model(pdb: &PDB, model_num: usize) -> PDB {
    let mut pdb_filtered = pdb.clone();
    if pdb_filtered.model_count() > 1 {
        // Find the model index (0-based) for the given model_num
        let model_idx = if model_num == 0 {
            0 // Use first model
        } else {
            pdb_filtered
                .models()
                .position(|m| m.serial_number() == model_num)
                .unwrap_or(0)
        };
        pdb_filtered.remove_models_except(&[model_idx]);
    }
    pdb_filtered
}

/// Common residue names for solvent molecules.
const SOLVENT_RESIDUES: &[&str] = &["HOH", "H2O", "D2O", "WAT", "TIP", "TIP3", "TIP4", "SPC"];

/// Common residue names for ions.
const ION_RESIDUES: &[&str] = &[
    "NA", "CL", "K", "CA", "MG", "ZN", "FE", "MN", "CU", "CO", "NI", "CD", "SO4", "PO4", "NO3",
    "ACE", "NH2",
];

/// Prepare a PDB structure for SASA calculation by removing solvent, ions, and hydrogens.
///
/// This function creates a clone of the PDB and removes:
/// - All hydrogen atoms
/// - Solvent molecules (HOH, H2O, WAT, etc.)
/// - Common ions (NA, CL, CA, MG, ZN, etc.)
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
///
/// # Returns
///
/// A filtered PDB structure suitable for SASA calculations.
pub(crate) fn prepare_pdb_for_sasa(pdb: &PDB) -> PDB {
    let mut pdb_prepared = pdb.clone();

    // Remove hydrogens, solvent, and ions
    pdb_prepared.remove_atoms_by(|atom| {
        // Remove hydrogens
        if atom.element() == Some(&Element::H) {
            return true;
        }
        false
    });

    // Remove entire residues that are solvent or ions
    // We need to check residue names using remove_residues_by
    pdb_prepared.remove_residues_by(|residue| {
        let resn = residue.name().unwrap_or("");
        SOLVENT_RESIDUES.contains(&resn) || ION_RESIDUES.contains(&resn)
    });

    pdb_prepared
}

/// Calculate solvent accessible surface area (SASA) for each atom in a PDB structure.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `num_threads` - Number of threads for parallel processing (1 for single-threaded, -1 for all cores)
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
/// let sasa_df = get_atom_sasa(&pdb, 1.4, 100, 0, 1);
/// println!("Calculated SASA for {} atoms", sasa_df.height());
/// ```
pub fn get_atom_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: isize,
) -> DataFrame {
    use rust_sasa::{Atom as SASAAtom, calculate_sasa_internal};

    // Prepare PDB: remove solvent, ions, and hydrogens
    let pdb_prepared = prepare_pdb_for_sasa(pdb);

    // If model_num is 0, we use the first model; otherwise use the specified model
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    // Calculate the SASA for each atom (excluding solvent, ions, hydrogens already removed)
    let atoms = pdb_filtered
        .atoms_with_hierarchy()
        .filter(|x| x.model().serial_number() == model_num)
        .map(|x| SASAAtom {
            position: [
                x.atom().pos().0 as f32,
                x.atom().pos().1 as f32,
                x.atom().pos().2 as f32,
            ],
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
    let atom_sasa = calculate_sasa_internal(&atoms, probe_radius, n_points, num_threads);

    // Create a DataFrame with the results
    let atom_annotations = pdb_filtered
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
    .sort(["atomi"], Default::default())
    .unwrap()
}

/// Calculate solvent accessible surface area (SASA) aggregated by residue.
///
/// Uses the rust-sasa SASAOptions API to compute SASA at the residue level.
/// Note there when there are multiple altlocs, only the first is considered.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `num_threads` - Number of threads for parallel processing (1 for single-threaded, -1 for all cores)
///
/// # Returns
///
/// A Polars DataFrame with columns:
/// - chain, resn, resi, insertion, sasa, is_polar
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_residue_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
/// let sasa_df = get_residue_sasa(&pdb, 1.4, 100, 0, 1);
/// println!("Calculated SASA for {} residues", sasa_df.height());
/// ```
pub fn get_residue_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: isize,
) -> DataFrame {
    use rust_sasa::{ResidueLevel, SASAOptions};

    // Prepare PDB: remove solvent, ions, and hydrogens, then filter by model
    let pdb_prepared = prepare_pdb_for_sasa(pdb);
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    let options = SASAOptions::<ResidueLevel>::new()
        .with_probe_radius(probe_radius)
        .with_n_points(n_points)
        .with_threads(num_threads)
        .with_allow_vdw_fallback(true);

    let result = options
        .process(&pdb_filtered)
        .expect("Failed to calculate residue-level SASA");

    df!(
        "chain" => result.iter().map(|r| r.chain_id.clone()).collect::<Vec<String>>(),
        "resn" => result.iter().map(|r| r.name.clone()).collect::<Vec<String>>(),
        "resi" => result.iter().map(|r| r.serial_number as i32).collect::<Vec<i32>>(),
        "insertion" => result.iter().map(|r| r.insertion_code.clone()).collect::<Vec<String>>(),
        "sasa" => result.iter().map(|r| r.value).collect::<Vec<f32>>(),
        "is_polar" => result.iter().map(|r| r.is_polar).collect::<Vec<bool>>(),
    )
    .unwrap()
    .sort(["chain", "resi", "insertion"], Default::default())
    .unwrap()
}

/// Calculate solvent accessible surface area (SASA) aggregated by chain.
///
/// Uses the rust-sasa SASAOptions API to compute SASA at the chain level.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `num_threads` - Number of threads for parallel processing (1 for single-threaded, -1 for all cores)
///
/// # Returns
///
/// A Polars DataFrame with columns:
/// - chain, sasa
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_chain_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
/// let sasa_df = get_chain_sasa(&pdb, 1.4, 100, 0, 1);
/// println!("Calculated SASA for {} chains", sasa_df.height());
/// ```
pub fn get_chain_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: isize,
) -> DataFrame {
    use rust_sasa::{ChainLevel, SASAOptions};

    // Prepare PDB: remove solvent, ions, and hydrogens, then filter by model
    let pdb_prepared = prepare_pdb_for_sasa(pdb);
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    let options = SASAOptions::<ChainLevel>::new()
        .with_probe_radius(probe_radius)
        .with_n_points(n_points)
        .with_threads(num_threads)
        .with_allow_vdw_fallback(true);

    let result = options
        .process(&pdb_filtered)
        .expect("Failed to calculate chain-level SASA");

    df!(
        "chain" => result.iter().map(|r| r.name.clone()).collect::<Vec<String>>(),
        "sasa" => result.iter().map(|r| r.value).collect::<Vec<f32>>(),
    )
    .unwrap()
    .sort(["chain"], Default::default())
    .unwrap()
}

/// Calculate the buried surface area (dSASA) at the interface between two chain groups.
///
/// The buried surface area is calculated as:
/// dSASA = (SASA_group1 + SASA_group2 - SASA_complex) / 2
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `groups` - Chain groups specification (e.g., "A,B/C,D")
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `num_threads` - Number of threads (-1 for all cores, 1 for single-threaded)
///
/// # Returns
///
/// The buried surface area at the interface in square Ångströms.
pub fn get_dsasa(
    pdb: &PDB,
    groups: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: isize,
) -> f32 {
    // Get all chains in the PDB
    let all_chains: HashSet<String> = pdb.chains().map(|c| c.id().to_string()).collect();

    // Parse groups
    let (group1_chains, group2_chains) = parse_groups(&all_chains, groups);

    // Get combined chains (union of both groups)
    let combined_group_chains: HashSet<String> =
        group1_chains.union(&group2_chains).cloned().collect();

    // Create PDB with only chains from both groups (remove unrelated chains)
    let mut pdb_combined = pdb.clone();
    pdb_combined.remove_chains_by(|chain| !combined_group_chains.contains(chain.id()));

    // Calculate SASA for the combined complex (only chains in groups)
    let combined_sasa = get_chain_sasa(
        &pdb_combined,
        probe_radius,
        n_points,
        model_num,
        num_threads,
    );

    let combined_total = sum_sasa(&combined_sasa);

    // Create PDB with only group1 chains and calculate SASA
    let mut pdb_group1 = pdb.clone();
    pdb_group1.remove_chains_by(|chain| !group1_chains.contains(chain.id()));

    let group1_sasa = get_chain_sasa(&pdb_group1, probe_radius, n_points, model_num, num_threads);

    let group1_total = sum_sasa(&group1_sasa);

    // Create PDB with only group2 chains and calculate SASA
    let mut pdb_group2 = pdb.clone();
    pdb_group2.remove_chains_by(|chain| !group2_chains.contains(chain.id()));

    let group2_sasa = get_chain_sasa(&pdb_group2, probe_radius, n_points, model_num, num_threads);

    let group2_total = sum_sasa(&group2_sasa);

    // Calculate buried surface area (dSASA)
    // dSASA = SASA_group1 + SASA_group2 - SASA_complex
    group1_total + group2_total - combined_total
}

/// Maximum solvent accessible surface area (MaxASA) values for amino acids.
///
/// Values are from Tien et al. (2013) "Maximum Allowed Solvent Accessibilities of Residues in Proteins"
/// PLOS ONE. These theoretical values represent the maximum possible SASA for each amino acid
/// in a Gly-X-Gly tripeptide.
///
/// Returns the MaxASA value in Å² for a given 3-letter amino acid code, or None if unknown.
pub fn get_max_asa(resn: &str) -> Option<f32> {
    match resn.to_uppercase().as_str() {
        "ALA" => Some(129.0),
        "ARG" => Some(274.0),
        "ASN" => Some(195.0),
        "ASP" => Some(193.0),
        "CYS" => Some(167.0),
        "GLU" => Some(223.0),
        "GLN" => Some(225.0),
        "GLY" => Some(104.0),
        "HIS" => Some(224.0),
        "ILE" => Some(197.0),
        "LEU" => Some(201.0),
        "LYS" => Some(236.0),
        "MET" => Some(224.0),
        "PHE" => Some(240.0),
        "PRO" => Some(159.0),
        "SER" => Some(155.0),
        "THR" => Some(172.0),
        "TRP" => Some(285.0),
        "TYR" => Some(263.0),
        "VAL" => Some(174.0),
        _ => None,
    }
}

/// Calculate relative solvent accessible surface area (RSA) for each residue.
///
/// RSA is calculated as the ratio of observed SASA to the maximum possible SASA
/// for each amino acid type, based on Tien et al. (2013) theoretical values.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `num_threads` - Number of threads (-1 for all cores, 1 for single-threaded)
///
/// # Returns
///
/// A Polars DataFrame with columns:
/// - chain, resn, resi, insertion, altloc, sasa, is_polar, max_sasa, relative_sasa
///
/// The relative_sasa column contains values between 0 and ~1 (can slightly exceed 1 due to
/// structural context), or null for non-standard amino acids.
pub fn get_relative_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: isize,
) -> DataFrame {
    // Get residue-level SASA
    let residue_sasa = get_residue_sasa(pdb, probe_radius, n_points, model_num, num_threads);

    // Calculate max_sasa and relative_sasa for each residue
    let resn_col = residue_sasa.column("resn").unwrap();
    let sasa_col = residue_sasa.column("sasa").unwrap();

    let max_sasa_values: Vec<Option<f32>> = resn_col
        .str()
        .unwrap()
        .into_iter()
        .map(|opt_resn| opt_resn.and_then(get_max_asa))
        .collect();

    let relative_sasa_values: Vec<Option<f32>> = sasa_col
        .f32()
        .unwrap()
        .into_iter()
        .zip(max_sasa_values.iter())
        .map(|(sasa_opt, max_opt)| match (sasa_opt, max_opt) {
            (Some(sasa), Some(max)) if *max > 0.0 => Some(sasa / max),
            _ => None,
        })
        .collect();

    // Add the relative_sasa column to the DataFrame
    let relative_sasa_series = Series::new("relative_sasa".into(), relative_sasa_values);

    residue_sasa
        .clone()
        .lazy()
        .with_columns([relative_sasa_series.lit()])
        .collect()
        .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::load_model;

    fn load_ubiquitin() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/1ubq.pdb");
        let (pdb, _) = load_model(&path);
        pdb
    }

    fn load_multi_chain() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/6bft.pdb");
        let (pdb, _) = load_model(&path);
        pdb
    }

    #[test]
    fn test_get_atom_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = get_atom_sasa(&pdb, 1.4, 100, 0, 1);

        // Check that we get results
        assert!(!df.is_empty(), "SASA DataFrame should not be empty");

        // Check that the expected columns exist
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            columns.contains(&"atomi".to_string()),
            "Should have 'atomi' column"
        );
        assert!(
            columns.contains(&"sasa".to_string()),
            "Should have 'sasa' column"
        );
        assert!(
            columns.contains(&"chain".to_string()),
            "Should have 'chain' column"
        );
        assert!(
            columns.contains(&"resn".to_string()),
            "Should have 'resn' column"
        );
        assert!(
            columns.contains(&"resi".to_string()),
            "Should have 'resi' column"
        );
        assert!(
            columns.contains(&"atomn".to_string()),
            "Should have 'atomn' column"
        );
    }

    #[test]
    fn test_get_atom_sasa_values_reasonable() {
        let pdb = load_ubiquitin();
        let df = get_atom_sasa(&pdb, 1.4, 100, 0, 1);

        // Get the SASA column and check values are non-negative
        let sasa_col = df.column("sasa").unwrap();
        let sasa_values: Vec<f32> = sasa_col.f32().unwrap().into_iter().flatten().collect();

        assert!(
            sasa_values.iter().all(|&v| v >= 0.0),
            "All SASA values should be non-negative"
        );

        // Check that some atoms have non-zero SASA (surface exposed)
        let non_zero_count = sasa_values.iter().filter(|&&v| v > 0.0).count();
        assert!(non_zero_count > 0, "Some atoms should have non-zero SASA");
    }

    #[test]
    fn test_get_residue_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = get_residue_sasa(&pdb, 1.4, 100, 0, 1);

        // Check that we get results
        assert!(!df.is_empty(), "Residue SASA DataFrame should not be empty");

        // Check that the expected columns exist
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            columns.contains(&"chain".to_string()),
            "Should have 'chain' column"
        );
        assert!(
            columns.contains(&"resn".to_string()),
            "Should have 'resn' column"
        );
        assert!(
            columns.contains(&"resi".to_string()),
            "Should have 'resi' column"
        );
        assert!(
            columns.contains(&"insertion".to_string()),
            "Should have 'insertion' column"
        );
        assert!(
            columns.contains(&"sasa".to_string()),
            "Should have 'sasa' column"
        );
        assert!(
            columns.contains(&"is_polar".to_string()),
            "Should have 'is_polar' column"
        );
    }

    #[test]
    fn test_get_residue_sasa_aggregation() {
        let pdb = load_ubiquitin();
        let atom_df = get_atom_sasa(&pdb, 1.4, 100, 0, 1);
        let residue_df = get_residue_sasa(&pdb, 1.4, 100, 0, 1);

        // There should be fewer rows in residue-level than atom-level
        assert!(
            residue_df.height() < atom_df.height(),
            "Residue-level should have fewer rows than atom-level: {} vs {}",
            residue_df.height(),
            atom_df.height()
        );

        // Total SASA at residue level should approximately match atom level
        // (may differ slightly due to different processing paths)
        let atom_total: f32 = sum_sasa(&atom_df);
        let residue_total: f32 = sum_sasa(&residue_df);

        // Allow for small differences due to potentially different filtering
        let ratio = residue_total / atom_total;
        assert!(
            ratio > 0.9 && ratio < 1.1,
            "Total SASA should be similar: atom={}, residue={}, ratio={}",
            atom_total,
            residue_total,
            ratio
        );
    }

    #[test]
    fn test_get_chain_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = get_chain_sasa(&pdb, 1.4, 100, 0, 1);

        // Check that we get results
        assert!(!df.is_empty(), "Chain SASA DataFrame should not be empty");

        // Check that the expected columns exist
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            columns.contains(&"chain".to_string()),
            "Should have 'chain' column"
        );
        assert!(
            columns.contains(&"sasa".to_string()),
            "Should have 'sasa' column"
        );
    }

    #[test]
    fn test_get_chain_sasa_single_chain() {
        let pdb = load_ubiquitin();
        let df = get_chain_sasa(&pdb, 1.4, 100, 0, 1);

        // Ubiquitin (1ubq) has a single chain A
        assert_eq!(df.height(), 1, "1ubq should have 1 chain");

        // Check that the chain is A
        let chain_col = df.column("chain").unwrap();
        let chain_id = chain_col.str().unwrap().get(0).unwrap();
        assert_eq!(chain_id, "A", "Chain should be A");
    }

    #[test]
    fn test_get_chain_sasa_multi_chain() {
        let pdb = load_multi_chain();
        let df = get_chain_sasa(&pdb, 1.4, 100, 0, 1);

        // 6bft should have multiple chains
        assert!(df.height() > 1, "6bft should have multiple chains");

        // Check that all SASA values are non-negative
        let sasa_col = df.column("sasa").unwrap();
        let sasa_values: Vec<f32> = sasa_col.f32().unwrap().into_iter().flatten().collect();

        assert!(
            sasa_values.iter().all(|&v| v >= 0.0),
            "All chain SASA values should be non-negative"
        );
    }

    #[test]
    fn test_sasa_probe_radius_effect() {
        let pdb = load_ubiquitin();

        // Smaller probe radius should result in larger SASA
        let small_probe = get_chain_sasa(&pdb, 1.0, 100, 0, 1);
        let large_probe = get_chain_sasa(&pdb, 2.0, 100, 0, 1);

        let small_sasa: f32 = small_probe
            .column("sasa")
            .unwrap()
            .f32()
            .unwrap()
            .get(0)
            .unwrap();
        let large_sasa: f32 = large_probe
            .column("sasa")
            .unwrap()
            .f32()
            .unwrap()
            .get(0)
            .unwrap();

        assert!(
            small_sasa > large_sasa,
            "Smaller probe radius should give larger SASA: {} vs {}",
            small_sasa,
            large_sasa
        );
    }

    #[test]
    fn test_sasa_regression_ubiquitin() {
        // Regression test to ensure SASA values remain consistent
        let pdb = load_ubiquitin();
        let df = get_chain_sasa(&pdb, 1.4, 100, 0, 1);

        let total_sasa: f32 = df.column("sasa").unwrap().f32().unwrap().get(0).unwrap();

        // Expected value from rust-sasa 0.9.0 with default parameters
        // This should be around 4813 Å² for ubiquitin
        let expected_sasa = 4813.0;
        let tolerance = 100.0; // Allow some tolerance for minor differences

        assert!(
            (total_sasa - expected_sasa).abs() < tolerance,
            "Ubiquitin total SASA should be around {} Å², got {} Å²",
            expected_sasa,
            total_sasa
        );
    }

    #[test]
    fn test_get_dsasa_returns_positive() {
        // 6bft has multiple chains with interfaces
        let pdb = load_multi_chain();

        // Calculate dSASA between groups A,B,C and G,H,L
        let dsasa = get_dsasa(&pdb, "A,B,C/G,H,L", 1.4, 100, 0, 1);

        // dSASA should be positive for an interface
        assert!(dsasa > 0.0, "dSASA should be positive, got {}", dsasa);
    }

    #[test]
    fn test_get_dsasa_interface_value() {
        // 6bft has multiple chains with interfaces
        let pdb = load_multi_chain();

        // Calculate dSASA between groups A,B,C and G,H,L
        let dsasa = get_dsasa(&pdb, "C/H,L", 1.4, 100, 0, 1);

        // Regression test: the dSASA should be around 1644-1665 Å²
        // as calculated from PyMOL and Rosetta InterfaceAnalyzer
        let expected_dsasa = 1650.0;
        let tolerance = 50.0; // Allow some tolerance

        assert!(
            (dsasa - expected_dsasa).abs() < tolerance,
            "6bft dSASA should be around {} Å², got {} Å²",
            expected_dsasa,
            dsasa
        );
    }

    #[test]
    fn test_get_dsasa_symmetric() {
        // dSASA should be the same regardless of which group is first
        let pdb = load_multi_chain();

        let dsasa1 = get_dsasa(&pdb, "A,B,C/G,H,L", 1.4, 100, 0, 1);
        let dsasa2 = get_dsasa(&pdb, "G,H,L/A,B,C", 1.4, 100, 0, 1);

        let diff = (dsasa1 - dsasa2).abs();
        assert!(
            diff < 1.0,
            "dSASA should be symmetric: {} vs {}",
            dsasa1,
            dsasa2
        );
    }

    #[test]
    fn test_get_relative_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = get_relative_sasa(&pdb, 1.4, 100, 0, 1);

        // Check that we get results
        assert!(
            !df.is_empty(),
            "Relative SASA DataFrame should not be empty"
        );

        // Check that the expected columns exist
        let columns: Vec<String> = df
            .get_column_names()
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert!(
            columns.contains(&"chain".to_string()),
            "Should have 'chain' column"
        );
        assert!(
            columns.contains(&"resn".to_string()),
            "Should have 'resn' column"
        );
        assert!(
            columns.contains(&"resi".to_string()),
            "Should have 'resi' column"
        );
        assert!(
            columns.contains(&"sasa".to_string()),
            "Should have 'sasa' column"
        );
        assert!(
            columns.contains(&"relative_sasa".to_string()),
            "Should have 'relative_sasa' column"
        );
    }

    #[test]
    fn test_get_relative_sasa_values_bounded() {
        let pdb = load_ubiquitin();
        let df = get_relative_sasa(&pdb, 1.4, 100, 0, 1);

        // Get relative_sasa values
        let rsa_values: Vec<f32> = df
            .column("relative_sasa")
            .unwrap()
            .f32()
            .unwrap()
            .into_iter()
            .flatten()
            .collect();

        // All values should be non-negative
        assert!(
            rsa_values.iter().all(|&v| v >= 0.0),
            "All relative_sasa values should be non-negative"
        );

        // Most values should be <= 1.0 (some may slightly exceed due to structural context)
        let below_threshold = rsa_values.iter().filter(|&&v| v <= 1.5).count();
        let ratio = below_threshold as f64 / rsa_values.len() as f64;
        assert!(
            ratio > 0.95,
            "Most relative_sasa values should be <= 1.5: ratio={}",
            ratio
        );
    }

    #[test]
    fn test_get_max_asa_standard_amino_acids() {
        // Test that all standard amino acids have MaxASA values
        let amino_acids = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ];

        for aa in amino_acids.iter() {
            let max_asa = get_max_asa(aa);
            assert!(max_asa.is_some(), "Should have MaxASA value for {}", aa);
            assert!(
                max_asa.unwrap() > 0.0,
                "MaxASA for {} should be positive",
                aa
            );
        }
    }

    #[test]
    fn test_get_max_asa_unknown_residue() {
        // Unknown residues should return None
        assert!(get_max_asa("XXX").is_none());
        assert!(get_max_asa("HOH").is_none());
        assert!(get_max_asa("").is_none());
    }
}
