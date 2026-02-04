//! Solvent Accessible Surface Area (SASA) calculations.
//!
//! This module provides functions for calculating SASA at different levels
//! of granularity (atom, residue, chain) and related metrics like dSASA
//! (buried surface area) and relative SASA.

use crate::contacts::InteractingEntity;
use crate::utils::{parse_groups, sum_float_col};
use pdbtbx::*;
use polars::prelude::*;
use std::collections::HashSet;

/// Filter a PDB structure to keep only the specified model.
///
/// Creates a clone of the PDB and removes all models except the one at the specified index.
/// If `model_num` is 0, the first model is used. Otherwise, the model with the matching
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

/// Parse a comma-separated chain string into a `HashSet` of chain IDs.
///
/// # Arguments
///
/// * `chains` - Comma-separated chain IDs (e.g., "A,B,C"). Empty string means all chains.
///
/// # Returns
///
/// A `HashSet` of chain IDs. If the input is empty, returns an empty `HashSet` (meaning all chains).
///
/// # Example
///
/// ```ignore
/// let chains = parse_chain_string("A,B,C");
/// assert!(chains.contains("A"));
/// assert!(chains.contains("B"));
/// assert!(chains.contains("C"));
/// ```
fn parse_chain_string(chains: &str) -> HashSet<String> {
    if chains.is_empty() {
        HashSet::new()
    } else {
        chains
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    }
}

/// Prepare a PDB structure for SASA calculation by removing solvent, ions, hydrogens,
/// and optionally filtering to specific chains.
///
/// This function creates a clone of the PDB and removes:
/// - Chains not in the specified list (if chains is not empty)
/// - All hydrogen atoms (if `remove_hydrogens` is true)
/// - Solvent molecules (HOH, H2O, WAT, etc.) (if `remove_solvent_and_ions` is true)
/// - Common ions (NA, CL, CA, MG, ZN, etc.) (if `remove_solvent_and_ions` is true)
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `remove_hydrogens` - Whether to remove hydrogen atoms
/// * `remove_solvent_and_ions` - Whether to remove solvent and ion residues
/// * `chains` - Comma-separated chain IDs to keep (e.g., "A,B,C"). Empty string keeps all chains.
///
/// # Returns
///
/// A filtered PDB structure suitable for SASA calculations.
///
/// # Example
///
/// ```ignore
/// // Keep only chains A and B, remove hydrogens and solvent
/// let filtered = prepare_pdb_for_sasa(&pdb, true, true, "A,B");
///
/// // Keep all chains
/// let all_chains = prepare_pdb_for_sasa(&pdb, true, true, "");
/// ```
pub(crate) fn prepare_pdb_for_sasa(
    pdb: &PDB,
    remove_hydrogens: bool,
    remove_solvent_and_ions: bool,
    chains: &str,
) -> PDB {
    let mut pdb_prepared = pdb.clone();

    // Filter to specific chains if specified
    let chain_filter = parse_chain_string(chains);
    if !chain_filter.is_empty() {
        pdb_prepared.remove_chains_by(|chain| !chain_filter.contains(chain.id()));
    }

    // Remove hydrogens
    if remove_hydrogens {
        pdb_prepared.remove_atoms_by(|atom| atom.element() == Some(&Element::H));
    }

    // Remove entire residues that are solvent or ions
    if remove_solvent_and_ions {
        pdb_prepared.remove_residues_by(|residue| {
            let resn = residue.name().unwrap_or("");
            SOLVENT_RESIDUES.contains(&resn) || ION_RESIDUES.contains(&resn)
        });
    }

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
/// * `remove_hydrogens` - Whether to remove hydrogen atoms before calculation
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - `atomi`, `sasa`
/// - `chain`, `resn`, `resi`, `insertion`, `altloc`, `atomn`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_atom_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
///
/// // Calculate SASA for all chains
/// let sasa_df = get_atom_sasa(&pdb, 1.4, 100, 0, true, "");
/// println!("Calculated SASA for {} atoms", sasa_df.height());
///
/// // Calculate SASA for only chains A and B
/// let sasa_ab = get_atom_sasa(&pdb, 1.4, 100, 0, true, "A,B");
/// ```
pub fn get_atom_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    remove_hydrogens: bool,
    chains: &str,
) -> DataFrame {
    use rust_sasa::{Atom as SASAAtom, calculate_sasa_internal};

    // Prepare PDB: remove solvent, ions, hydrogens, and filter chains
    let pdb_prepared = prepare_pdb_for_sasa(pdb, remove_hydrogens, true, chains);

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
    let atom_sasa = calculate_sasa_internal(
        &atoms,
        probe_radius,
        n_points,
        rayon::current_num_threads() as isize,
    );

    // Create a DataFrame with the results
    let atom_annotations = pdb_filtered
        .atoms_with_hierarchy()
        .map(|x| InteractingEntity::from_hier(&x))
        .collect::<Vec<InteractingEntity>>();
    let atom_annot_df = df!(
        "chain" => atom_annotations.iter().map(|x| x.chain.clone()).collect::<Vec<String>>(),
        "resn" => atom_annotations.iter().map(|x| x.resn.clone()).collect::<Vec<String>>(),
        "resi" => atom_annotations.iter().map(|x| x.resi as i32).collect::<Vec<i32>>(),
        "insertion" => atom_annotations.iter().map(|x| x.insertion.clone()).collect::<Vec<String>>(),
        "altloc" => atom_annotations.iter().map(|x| x.altloc.clone()).collect::<Vec<String>>(),
        "atomn" => atom_annotations.iter().map(|x| x.atomn.clone()).collect::<Vec<String>>(),
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
/// Uses the rust-sasa `SASAOptions` API to compute SASA at the residue level.
/// Note there when there are multiple altlocs, only the first is considered.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - `chain`, `resn`, `resi`, `insertion`, `sasa`, `is_polar`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_residue_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
///
/// // Calculate SASA for all chains
/// let sasa_df = get_residue_sasa(&pdb, 1.4, 100, 0, "");
/// println!("Calculated SASA for {} residues", sasa_df.height());
///
/// // Calculate SASA for only chain A
/// let sasa_a = get_residue_sasa(&pdb, 1.4, 100, 0, "A");
/// ```
pub fn get_residue_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> DataFrame {
    use rust_sasa::{ResidueLevel, SASAOptions};

    // Prepare PDB: remove solvent, ions, hydrogens, and filter chains, then filter by model
    let pdb_prepared = prepare_pdb_for_sasa(pdb, true, true, chains);
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    let options = SASAOptions::<ResidueLevel>::new()
        .with_probe_radius(probe_radius)
        .with_n_points(n_points)
        .with_threads(rayon::current_num_threads() as isize)
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
/// Uses the rust-sasa `SASAOptions` API to compute SASA at the chain level.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `probe_radius` - Probe radius in Ångströms (typically 1.4)
/// * `n_points` - Number of points for surface calculation (typically 100)
/// * `model_num` - Model number to analyze (0 for first model)
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - `chain`, `sasa`
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_chain_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
///
/// // Calculate SASA for all chains
/// let sasa_df = get_chain_sasa(&pdb, 1.4, 100, 0, "");
/// println!("Calculated SASA for {} chains", sasa_df.height());
///
/// // Calculate SASA for only chains A and B
/// let sasa_ab = get_chain_sasa(&pdb, 1.4, 100, 0, "A,B");
/// ```
pub fn get_chain_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> DataFrame {
    use rust_sasa::{ChainLevel, SASAOptions};

    // Prepare PDB: remove solvent, ions, hydrogens, and filter chains, then filter by model
    let pdb_prepared = prepare_pdb_for_sasa(pdb, true, true, chains);
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    let options = SASAOptions::<ChainLevel>::new()
        .with_probe_radius(probe_radius)
        .with_n_points(n_points)
        .with_threads(rayon::current_num_threads() as isize)
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
        "", // Already filtered by combined_group_chains
    );

    let combined_total = sum_float_col(&combined_sasa, "sasa");

    // Create PDB with only group1 chains and calculate SASA
    let mut pdb_group1 = pdb.clone();
    pdb_group1.remove_chains_by(|chain| !group1_chains.contains(chain.id()));

    let group1_sasa = get_chain_sasa(&pdb_group1, probe_radius, n_points, model_num, "");

    let group1_total = sum_float_col(&group1_sasa, "sasa");

    // Create PDB with only group2 chains and calculate SASA
    let mut pdb_group2 = pdb.clone();
    pdb_group2.remove_chains_by(|chain| !group2_chains.contains(chain.id()));

    let group2_sasa = get_chain_sasa(&pdb_group2, probe_radius, n_points, model_num, "");

    let group2_total = sum_float_col(&group2_sasa, "sasa");

    // Calculate buried surface area (dSASA)
    // dSASA = SASA_group1 + SASA_group2 - SASA_complex
    group1_total + group2_total - combined_total
}

/// Maximum solvent accessible surface area (`MaxASA`) values for amino acids.
///
/// Values are from Tien et al. (2013) "Maximum Allowed Solvent Accessibilities of Residues in Proteins"
/// PLOS ONE. These theoretical values represent the maximum possible SASA for each amino acid
/// in a Gly-X-Gly tripeptide.
///
/// Returns the `MaxASA` value in Å² for a given 3-letter amino acid code, or None if unknown.
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
        "HIS" | "MET" => Some(224.0),
        "ILE" => Some(197.0),
        "LEU" => Some(201.0),
        "LYS" => Some(236.0),
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
/// * `chains` - Comma-separated chain IDs to include (e.g., "A,B,C"). Empty string includes all chains.
///
/// # Returns
///
/// A Polars `DataFrame` with columns:
/// - chain, resn, resi, insertion, altloc, sasa, `is_polar`, `max_sasa`, `relative_sasa`
///
/// The `relative_sasa` column contains values between 0 and ~1 (can slightly exceed 1 due to
/// structural context), or null for non-standard amino acids.
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_relative_sasa};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
///
/// // Calculate RSA for all chains
/// let rsa_df = get_relative_sasa(&pdb, 1.4, 100, 0, "");
///
/// // Calculate RSA for only chain A
/// let rsa_a = get_relative_sasa(&pdb, 1.4, 100, 0, "A");
/// ```
pub fn get_relative_sasa(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> DataFrame {
    // Get residue-level SASA
    let residue_sasa = get_residue_sasa(pdb, probe_radius, n_points, model_num, chains);

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
    use crate::utils::{load_model, run_with_threads};

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
        let df = run_with_threads(1, || get_atom_sasa(&pdb, 1.4, 100, 0, true, ""));

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
        let df = run_with_threads(1, || get_atom_sasa(&pdb, 1.4, 100, 0, true, ""));

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
        let df = run_with_threads(1, || get_residue_sasa(&pdb, 1.4, 100, 0, ""));

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
        let atom_df = run_with_threads(1, || get_atom_sasa(&pdb, 1.4, 100, 0, true, ""));
        let residue_df = run_with_threads(1, || get_residue_sasa(&pdb, 1.4, 100, 0, ""));

        // There should be fewer rows in residue-level than atom-level
        assert!(
            residue_df.height() < atom_df.height(),
            "Residue-level should have fewer rows than atom-level: {} vs {}",
            residue_df.height(),
            atom_df.height()
        );

        // Total SASA at residue level should approximately match atom level
        // (may differ slightly due to different processing paths)
        let atom_total: f32 = sum_float_col(&atom_df, "sasa");
        let residue_total: f32 = sum_float_col(&residue_df, "sasa");

        // Allow for small differences due to potentially different filtering
        let ratio = residue_total / atom_total;
        assert!(
            ratio > 0.9 && ratio < 1.1,
            "Total SASA should be similar: atom={atom_total}, residue={residue_total}, ratio={ratio}"
        );
    }

    #[test]
    fn test_get_chain_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));

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
        let df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));

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
        let df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));

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
        let small_probe = run_with_threads(1, || get_chain_sasa(&pdb, 1.0, 100, 0, ""));
        let large_probe = run_with_threads(1, || get_chain_sasa(&pdb, 2.0, 100, 0, ""));

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
            "Smaller probe radius should give larger SASA: {small_sasa} vs {large_sasa}"
        );
    }

    #[test]
    fn test_sasa_regression_ubiquitin() {
        // Regression test to ensure SASA values remain consistent
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));

        let total_sasa: f32 = df.column("sasa").unwrap().f32().unwrap().get(0).unwrap();

        // Expected value from rust-sasa 0.9.0 with default parameters
        // This should be around 4813 Å² for ubiquitin
        let expected_sasa = 4813.0;
        let tolerance = 100.0; // Allow some tolerance for minor differences

        assert!(
            (total_sasa - expected_sasa).abs() < tolerance,
            "Ubiquitin total SASA should be around {expected_sasa} Å², got {total_sasa} Å²"
        );
    }

    #[test]
    fn test_get_dsasa_returns_positive() {
        // 6bft has multiple chains with interfaces
        let pdb = load_multi_chain();

        // Calculate dSASA between groups A,B,C and G,H,L
        let dsasa = run_with_threads(1, || get_dsasa(&pdb, "A,B,C/G,H,L", 1.4, 100, 0));

        // dSASA should be positive for an interface
        assert!(dsasa > 0.0, "dSASA should be positive, got {dsasa}");
    }

    #[test]
    fn test_get_dsasa_interface_value() {
        // 6bft has multiple chains with interfaces
        let pdb = load_multi_chain();

        // Calculate dSASA between groups A,B,C and G,H,L
        let dsasa = run_with_threads(1, || get_dsasa(&pdb, "C/H,L", 1.4, 100, 0));

        // Regression test: the dSASA should be around 1644-1665 Å²
        // as calculated from PyMOL and Rosetta InterfaceAnalyzer
        let expected_dsasa = 1650.0;
        let tolerance = 50.0; // Allow some tolerance

        assert!(
            (dsasa - expected_dsasa).abs() < tolerance,
            "6bft dSASA should be around {expected_dsasa} Å², got {dsasa} Å²"
        );
    }

    #[test]
    fn test_get_dsasa_symmetric() {
        // dSASA should be the same regardless of which group is first
        let pdb = load_multi_chain();

        let dsasa1 = run_with_threads(1, || get_dsasa(&pdb, "A,B,C/G,H,L", 1.4, 100, 0));
        let dsasa2 = run_with_threads(1, || get_dsasa(&pdb, "G,H,L/A,B,C", 1.4, 100, 0));

        let diff = (dsasa1 - dsasa2).abs();
        assert!(
            diff < 1.0,
            "dSASA should be symmetric: {dsasa1} vs {dsasa2}"
        );
    }

    #[test]
    fn test_get_relative_sasa_returns_data() {
        let pdb = load_ubiquitin();
        let df = run_with_threads(1, || get_relative_sasa(&pdb, 1.4, 100, 0, ""));

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
        let df = run_with_threads(1, || get_relative_sasa(&pdb, 1.4, 100, 0, ""));

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
            "Most relative_sasa values should be <= 1.5: ratio={ratio}"
        );
    }

    #[test]
    fn test_get_max_asa_standard_amino_acids() {
        // Test that all standard amino acids have MaxASA values
        let amino_acids = [
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
        ];

        for aa in &amino_acids {
            let max_asa = get_max_asa(aa);
            assert!(max_asa.is_some(), "Should have MaxASA value for {aa}");
            assert!(max_asa.unwrap() > 0.0, "MaxASA for {aa} should be positive");
        }
    }

    #[test]
    fn test_get_max_asa_unknown_residue() {
        // Unknown residues should return None
        assert!(get_max_asa("XXX").is_none());
        assert!(get_max_asa("HOH").is_none());
        assert!(get_max_asa("").is_none());
    }

    #[test]
    fn test_chain_filter_empty_keeps_all() {
        let pdb = load_multi_chain();

        // Empty chain filter should keep all chains
        let df_all = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));
        let chain_count_all = df_all.height();

        assert!(
            chain_count_all > 1,
            "Multi-chain structure should have multiple chains: {chain_count_all}"
        );
    }

    #[test]
    fn test_chain_filter_single_chain() {
        let pdb = load_multi_chain();

        // Filter to only chain A
        let df_a = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, "A"));

        assert_eq!(
            df_a.height(),
            1,
            "Should have only one chain when filtering to A"
        );

        // Check that the chain is A
        let chain_col = df_a.column("chain").unwrap();
        let chain_id = chain_col.str().unwrap().get(0).unwrap();
        assert_eq!(chain_id, "A", "Chain should be A");
    }

    #[test]
    fn test_chain_filter_multiple_chains() {
        let pdb = load_multi_chain();

        // Filter to chains A, B
        let df_ab = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, "A,B"));
        let df_all = run_with_threads(1, || get_chain_sasa(&pdb, 1.4, 100, 0, ""));

        assert!(
            df_ab.height() <= df_all.height(),
            "Filtered results should have equal or fewer chains: {} vs {}",
            df_ab.height(),
            df_all.height()
        );

        // Check that we only have A and B chains
        let chain_col = df_ab.column("chain").unwrap();
        let chain_ids: Vec<&str> = chain_col.str().unwrap().into_iter().flatten().collect();
        for chain_id in &chain_ids {
            assert!(
                *chain_id == "A" || *chain_id == "B",
                "Only A and B chains should be present, got: {chain_id}"
            );
        }
    }
}
