//! Solvent Accessible Surface Area (SASA) calculations.
//!
//! This module provides functions for calculating SASA at different levels
//! of granularity (atom, residue, chain) and related metrics like dSASA
//! (buried surface area) and relative SASA.

use crate::contacts::InteractingEntity;
use crate::structure::{parse_groups, prepare_structure, prepare_structure_with_chains};
use pdbtbx::*;
use polars::prelude::*;
use std::collections::HashSet;

#[derive(Clone, Debug)]
pub(crate) struct AtomSasaRecord {
    pub(crate) chain: String,
    pub(crate) resn: String,
    pub(crate) resi: i32,
    pub(crate) insertion: String,
    pub(crate) altloc: String,
    pub(crate) atomn: String,
    pub(crate) atomi: i32,
    pub(crate) sasa: f32,
}

#[derive(Clone, Debug)]
pub(crate) struct ResidueSasaRecord {
    pub(crate) chain: String,
    pub(crate) resn: String,
    pub(crate) resi: i32,
    pub(crate) insertion: String,
    pub(crate) sasa: f32,
    pub(crate) is_polar: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct ChainSasaRecord {
    pub(crate) chain: String,
    pub(crate) sasa: f32,
}

#[derive(Clone, Debug)]
pub(crate) struct RelativeSasaRecord {
    pub(crate) chain: String,
    pub(crate) resn: String,
    pub(crate) resi: i32,
    pub(crate) insertion: String,
    pub(crate) sasa: f32,
    pub(crate) is_polar: bool,
    pub(crate) relative_sasa: Option<f32>,
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
    atom_sasa_records_to_dataframe(&calculate_atom_sasa_records(
        pdb,
        probe_radius,
        n_points,
        model_num,
        remove_hydrogens,
        chains,
    ))
}

pub(crate) fn calculate_atom_sasa_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    remove_hydrogens: bool,
    chains: &str,
) -> Vec<AtomSasaRecord> {
    use rust_sasa::{Atom as SASAAtom, calculate_sasa_internal};

    let pdb_filtered = prepare_structure(pdb, model_num, remove_hydrogens, chains);

    let atom_hierarchy = pdb_filtered.atoms_with_hierarchy().collect::<Vec<_>>();
    let atoms = atom_hierarchy
        .iter()
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

    atom_hierarchy
        .iter()
        .zip(atom_sasa)
        .map(|(hier, sasa)| {
            let entity = InteractingEntity::from_hier(hier);
            AtomSasaRecord {
                chain: entity.chain,
                resn: entity.resn,
                resi: entity.resi as i32,
                insertion: entity.insertion,
                altloc: entity.altloc,
                atomn: entity.atomn,
                atomi: entity.atomi as i32,
                sasa,
            }
        })
        .collect()
}

pub(crate) fn atom_sasa_records_to_dataframe(records: &[AtomSasaRecord]) -> DataFrame {
    df!(
        "atomi" => records.iter().map(|x| x.atomi).collect::<Vec<i32>>(),
        "sasa" => records.iter().map(|x| x.sasa).collect::<Vec<f32>>(),
        "chain" => records.iter().map(|x| x.chain.clone()).collect::<Vec<String>>(),
        "resn" => records.iter().map(|x| x.resn.clone()).collect::<Vec<String>>(),
        "resi" => records.iter().map(|x| x.resi).collect::<Vec<i32>>(),
        "insertion" => records.iter().map(|x| x.insertion.clone()).collect::<Vec<String>>(),
        "altloc" => records.iter().map(|x| x.altloc.clone()).collect::<Vec<String>>(),
        "atomn" => records.iter().map(|x| x.atomn.clone()).collect::<Vec<String>>(),
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
    residue_sasa_records_to_dataframe(&calculate_residue_sasa_records(
        pdb,
        probe_radius,
        n_points,
        model_num,
        chains,
    ))
}

pub(crate) fn calculate_residue_sasa_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> Vec<ResidueSasaRecord> {
    use rust_sasa::{ResidueLevel, SASAOptions};

    let pdb_filtered = prepare_structure(pdb, model_num, true, chains);

    let options = SASAOptions::<ResidueLevel>::new()
        .with_probe_radius(probe_radius)
        .with_n_points(n_points)
        .with_threads(rayon::current_num_threads() as isize)
        .with_allow_vdw_fallback(true);

    let result = options
        .process(&pdb_filtered)
        .expect("Failed to calculate residue-level SASA");

    let mut records = result
        .iter()
        .map(|r| ResidueSasaRecord {
            chain: r.chain_id.clone(),
            resn: r.name.clone(),
            resi: r.serial_number as i32,
            insertion: r.insertion_code.clone(),
            sasa: r.value,
            is_polar: r.is_polar,
        })
        .collect::<Vec<_>>();
    records.sort_by(|a, b| (&a.chain, a.resi, &a.insertion).cmp(&(&b.chain, b.resi, &b.insertion)));
    records
}

pub(crate) fn residue_sasa_records_to_dataframe(records: &[ResidueSasaRecord]) -> DataFrame {
    df!(
        "chain" => records.iter().map(|r| r.chain.clone()).collect::<Vec<String>>(),
        "resn" => records.iter().map(|r| r.resn.clone()).collect::<Vec<String>>(),
        "resi" => records.iter().map(|r| r.resi).collect::<Vec<i32>>(),
        "insertion" => records.iter().map(|r| r.insertion.clone()).collect::<Vec<String>>(),
        "sasa" => records.iter().map(|r| r.sasa).collect::<Vec<f32>>(),
        "is_polar" => records.iter().map(|r| r.is_polar).collect::<Vec<bool>>(),
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
    chain_sasa_records_to_dataframe(&calculate_chain_sasa_records(
        pdb,
        probe_radius,
        n_points,
        model_num,
        chains,
    ))
}

pub(crate) fn calculate_chain_sasa_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> Vec<ChainSasaRecord> {
    use rust_sasa::{ChainLevel, SASAOptions};

    let pdb_filtered = prepare_structure(pdb, model_num, true, chains);

    let options = SASAOptions::<ChainLevel>::new()
        .with_probe_radius(probe_radius)
        .with_n_points(n_points)
        .with_threads(rayon::current_num_threads() as isize)
        .with_allow_vdw_fallback(true);

    let result = options
        .process(&pdb_filtered)
        .expect("Failed to calculate chain-level SASA");

    let mut records = result
        .iter()
        .map(|r| ChainSasaRecord {
            chain: r.name.clone(),
            sasa: r.value,
        })
        .collect::<Vec<_>>();
    records.sort_by(|a, b| a.chain.cmp(&b.chain));
    records
}

fn calculate_chain_sasa_records_for_chain_set(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &HashSet<String>,
) -> Vec<ChainSasaRecord> {
    use rust_sasa::{ChainLevel, SASAOptions};

    let pdb_filtered = prepare_structure_with_chains(pdb, model_num, true, chains);

    let options = SASAOptions::<ChainLevel>::new()
        .with_probe_radius(probe_radius)
        .with_n_points(n_points)
        .with_threads(rayon::current_num_threads() as isize)
        .with_allow_vdw_fallback(true);

    let result = options
        .process(&pdb_filtered)
        .expect("Failed to calculate chain-level SASA");

    let mut records = result
        .iter()
        .map(|r| ChainSasaRecord {
            chain: r.name.clone(),
            sasa: r.value,
        })
        .collect::<Vec<_>>();
    records.sort_by(|a, b| a.chain.cmp(&b.chain));
    records
}

pub(crate) fn chain_sasa_records_to_dataframe(records: &[ChainSasaRecord]) -> DataFrame {
    df!(
        "chain" => records.iter().map(|r| r.chain.clone()).collect::<Vec<String>>(),
        "sasa" => records.iter().map(|r| r.sasa).collect::<Vec<f32>>(),
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

    let combined_sasa = calculate_chain_sasa_records_for_chain_set(
        pdb,
        probe_radius,
        n_points,
        model_num,
        &combined_group_chains,
    );
    let combined_total = sum_chain_sasa(&combined_sasa);

    let group1_sasa = calculate_chain_sasa_records_for_chain_set(
        pdb,
        probe_radius,
        n_points,
        model_num,
        &group1_chains,
    );
    let group1_total = sum_chain_sasa(&group1_sasa);

    let group2_sasa = calculate_chain_sasa_records_for_chain_set(
        pdb,
        probe_radius,
        n_points,
        model_num,
        &group2_chains,
    );
    let group2_total = sum_chain_sasa(&group2_sasa);

    group1_total + group2_total - combined_total
}

fn sum_chain_sasa(records: &[ChainSasaRecord]) -> f32 {
    records.iter().map(|record| record.sasa).sum()
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
    relative_sasa_records_to_dataframe(&calculate_relative_sasa_records(
        pdb,
        probe_radius,
        n_points,
        model_num,
        chains,
    ))
}

fn calculate_relative_sasa_records(
    pdb: &PDB,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    chains: &str,
) -> Vec<RelativeSasaRecord> {
    calculate_residue_sasa_records(pdb, probe_radius, n_points, model_num, chains)
        .into_iter()
        .map(|record| {
            let relative_sasa = get_max_asa(&record.resn)
                .filter(|max_sasa| *max_sasa > 0.0)
                .map(|max_sasa| record.sasa / max_sasa);

            RelativeSasaRecord {
                chain: record.chain,
                resn: record.resn,
                resi: record.resi,
                insertion: record.insertion,
                sasa: record.sasa,
                is_polar: record.is_polar,
                relative_sasa,
            }
        })
        .collect()
}

fn relative_sasa_records_to_dataframe(records: &[RelativeSasaRecord]) -> DataFrame {
    df!(
        "chain" => records.iter().map(|r| r.chain.clone()).collect::<Vec<String>>(),
        "resn" => records.iter().map(|r| r.resn.clone()).collect::<Vec<String>>(),
        "resi" => records.iter().map(|r| r.resi).collect::<Vec<i32>>(),
        "insertion" => records.iter().map(|r| r.insertion.clone()).collect::<Vec<String>>(),
        "sasa" => records.iter().map(|r| r.sasa).collect::<Vec<f32>>(),
        "is_polar" => records.iter().map(|r| r.is_polar).collect::<Vec<bool>>(),
        "relative_sasa" => records.iter().map(|r| r.relative_sasa).collect::<Vec<Option<f32>>>(),
    )
    .unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{load_model, run_with_threads, sum_float_col};

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
        assert!(df.height() > 0, "SASA DataFrame should not be empty");

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
        let sasa_values: Vec<f32> = sasa_col.f32().unwrap().iter().flatten().collect();

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
        assert!(
            df.height() > 0,
            "Residue SASA DataFrame should not be empty"
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
        assert!(df.height() > 0, "Chain SASA DataFrame should not be empty");

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
        let sasa_values: Vec<f32> = sasa_col.f32().unwrap().iter().flatten().collect();

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
            df.height() > 0,
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
            .iter()
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
        let chain_ids: Vec<&str> = chain_col.str().unwrap().iter().flatten().collect();
        for chain_id in &chain_ids {
            assert!(
                *chain_id == "A" || *chain_id == "B",
                "Only A and B chains should be present, got: {chain_id}"
            );
        }
    }
}
