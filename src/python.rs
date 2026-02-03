//! Python bindings for the arpeggia library using PyO3.
//!
//! This module provides Python-friendly wrappers around the core Rust functions,
//! converting Polars DataFrames to Python using pyo3-polars for efficient zero-copy data transfer.

use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

/// Convert thread count from Python usize to rust-sasa isize.
fn get_num_threads(num_threads: usize) -> isize {
    crate::get_num_threads(num_threads)
}

/// Load a PDB or mmCIF file and calculate atomic and ring contacts.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     groups (str, optional): Chain groups specification. Defaults to "/" (all-to-all).
///         Examples: "A,B/C,D" for chains A,B vs C,D; "A/" for chain A vs all others.
///     vdw_comp (float, optional): VdW radii compensation factor. Defaults to 0.1.
///     dist_cutoff (float, optional): Distance cutoff for neighbor searches in Ångströms. Defaults to 6.5.
///     ignore_zero_occupancy (bool, optional): If True, ignore atoms with zero occupancy. Defaults to False.
///
/// Returns:
///     polars.DataFrame: A DataFrame containing all identified contacts with columns:
///         - model, interaction, distance
///         - from_chain, from_resn, from_resi, from_insertion, from_altloc, from_atomn, from_atomi
///         - to_chain, to_resn, to_resi, to_insertion, to_altloc, to_atomn, to_atomi
///         - sc_centroid_dist, sc_dihedral, sc_centroid_angle
///
/// Example:
///     >>> import arpeggia
///     >>> contacts = arpeggia.contacts("structure.pdb", groups="/", vdw_comp=0.1)
///     >>> print(f"Found {len(contacts)} contacts")
#[pyfunction]
#[pyo3(signature = (input_file, groups="/", vdw_comp=0.1, dist_cutoff=6.5, ignore_zero_occupancy=false))]
fn contacts(
    input_file: String,
    groups: &str,
    vdw_comp: f64,
    dist_cutoff: f64,
    ignore_zero_occupancy: bool,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let (mut pdb, _warnings) = crate::load_model(&input_file);

    // Filter out atoms with zero occupancy if requested
    if ignore_zero_occupancy {
        pdb.remove_atoms_by(|atom| atom.occupancy() == 0.0);
    }

    // Get contacts
    let df = crate::get_contacts(&pdb, groups, vdw_comp, dist_cutoff);

    // Convert to PyDataFrame for Python
    Ok(PyDataFrame(df))
}

/// Load a PDB or mmCIF file and calculate solvent accessible surface area (SASA).
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     level (str, optional): Aggregation level for SASA calculation. Options:
///         - "atom": Calculate SASA for each atom (default)
///         - "residue": Aggregate SASA by residue
///         - "chain": Aggregate SASA by chain
///     probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
///     n_points (int, optional): Number of points for surface calculation. Defaults to 100.
///     model_num (int, optional): Model number to analyze (0 for first model). Defaults to 0.
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///     chains (str, optional): Comma-separated chain IDs to include (e.g., "A,B,C").
///         If empty, includes all chains. Defaults to "".
///
/// Returns:
///     polars.DataFrame: A DataFrame with SASA values. Columns depend on the level:
///         - atom: atomi, sasa, chain, resn, resi, insertion, altloc, atomn
///         - residue: chain, resn, resi, insertion, altloc, sasa, is_polar
///         - chain: chain, sasa
///
/// Example:
///     >>> import arpeggia
///     >>> # Atom-level SASA for all chains
///     >>> sasa_df = arpeggia.sasa("structure.pdb", level="atom")
///     >>> print(f"Calculated SASA for {len(sasa_df)} atoms")
///     >>>
///     >>> # Residue-level SASA for only chains A and B
///     >>> residue_sasa = arpeggia.sasa("structure.pdb", level="residue", chains="A,B")
///     >>> print(f"Calculated SASA for {len(residue_sasa)} residues")
///     >>>
///     >>> # Chain-level SASA
///     >>> chain_sasa = arpeggia.sasa("structure.pdb", level="chain")
///     >>> print(f"Calculated SASA for {len(chain_sasa)} chains")
#[pyfunction]
#[pyo3(signature = (input_file, level="atom", probe_radius=1.4, n_points=100, model_num=0, num_threads=1, chains=""))]
fn sasa(
    input_file: String,
    level: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: usize,
    chains: &str,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Convert num_threads for rust-sasa
    let threads = get_num_threads(num_threads);

    // Get SASA based on level
    let df = match level.to_lowercase().as_str() {
        "atom" => crate::get_atom_sasa(&pdb, probe_radius, n_points, model_num, threads, true, chains),
        "residue" => crate::get_residue_sasa(&pdb, probe_radius, n_points, model_num, threads, chains),
        "chain" => crate::get_chain_sasa(&pdb, probe_radius, n_points, model_num, threads, chains),
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid level '{}'. Must be one of: 'atom', 'residue', 'chain'",
                level
            )));
        }
    };

    // Convert to PyDataFrame for Python
    Ok(PyDataFrame(df))
}

/// Load a PDB or mmCIF file and calculate buried surface area at the interface between chain groups.
///
/// The buried surface area (dSASA) is calculated as:
/// dSASA = SASA_group1 + SASA_group2 - SASA_complex
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     groups (str): Chain groups specification for interface calculation.
///         Format: "A,B/C,D" where chains A,B form one side and C,D form the other side.
///     probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
///     n_points (int, optional): Number of points for surface calculation. Defaults to 100.
///     model_num (int, optional): Model number to analyze (0 for first model). Defaults to 0.
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///
/// Returns:
///     float: The buried surface area at the interface in square Ångströms.
///
/// Example:
///     >>> import arpeggia
///     >>> bsa = arpeggia.dsasa("structure.pdb", groups="A,B/C,D")
///     >>> print(f"Buried surface area: {bsa:.2f} Å²")
#[pyfunction]
#[pyo3(signature = (input_file, groups, probe_radius=1.4, n_points=100, model_num=0, num_threads=1))]
fn dsasa(
    input_file: String,
    groups: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: usize,
) -> PyResult<f32> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Convert num_threads for rust-sasa
    let threads = get_num_threads(num_threads);

    // Use the library function to calculate dSASA
    let result = crate::get_dsasa(&pdb, groups, probe_radius, n_points, model_num, threads);

    if result <= 0.0 {
        // dSASA can be 0 for non-interacting groups, but we should warn about potential issues
        // A negative dSASA would indicate an error in the calculation
        if result < 0.0 {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "Negative dSASA calculated. Please check the input file and chain groups.",
            ));
        }
    }

    Ok(result)
}

/// Load a PDB or mmCIF file and extract sequences for all chains.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///
/// Returns:
///     dict[str, str]: A dictionary mapping chain IDs to their sequences
///
/// Example:
///     >>> import arpeggia
///     >>> sequences = arpeggia.sequences("structure.pdb")
///     >>> for chain_id, seq in sequences.items():
///     ...     print(f"Chain {chain_id}: {seq}")
#[pyfunction]
fn pdb2seq(input_file: String) -> PyResult<std::collections::HashMap<String, String>> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Get sequences
    let seqs = crate::get_sequences(&pdb);

    Ok(seqs)
}

/// Load a PDB or mmCIF file and calculate relative SASA (RSA) for each residue.
///
/// RSA is calculated as the ratio of observed SASA to the maximum possible SASA
/// for each amino acid type, based on Tien et al. (2013) theoretical values.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
///     n_points (int, optional): Number of points for surface calculation. Defaults to 100.
///     model_num (int, optional): Model number to analyze (0 for first model). Defaults to 0.
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///     chains (str, optional): Comma-separated chain IDs to include (e.g., "A,B,C").
///         If empty, includes all chains. Defaults to "".
///
/// Returns:
///     polars.DataFrame: A DataFrame with relative SASA values for each residue with columns:
///         - chain, resn, resi, insertion, altloc, sasa, is_polar, relative_sasa
///
/// Example:
///     >>> import arpeggia
///     >>> # RSA for all chains
///     >>> rsa = arpeggia.relative_sasa("structure.pdb", probe_radius=1.4)
///     >>> print(f"Calculated RSA for {len(rsa)} residues")
///     >>>
///     >>> # RSA for only chain A
///     >>> rsa_a = arpeggia.relative_sasa("structure.pdb", chains="A")
#[pyfunction]
#[pyo3(signature = (input_file, probe_radius=1.4, n_points=100, model_num=0, num_threads=1, chains=""))]
fn relative_sasa(
    input_file: String,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    num_threads: usize,
    chains: &str,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Convert num_threads for rust-sasa
    let threads = get_num_threads(num_threads);

    // Get relative SASA
    let df = crate::get_relative_sasa(&pdb, probe_radius, n_points, model_num, threads, chains);

    // Convert to PyDataFrame for Python
    Ok(PyDataFrame(df))
}

/// Load a PDB or mmCIF file and calculate Spatial Aggregation Propensity (SAP) scores.
///
/// The SAP score quantifies the aggregation propensity by combining the solvent-accessible
/// hydrophobic surface area of neighboring residues. It was developed by Chennamsetty et al.
/// and is described in "Developability Index: A Rapid In Silico Tool for the Screening of
/// Antibody Aggregation Propensity" (J Pharm Sci, 2012).
///
/// The formula is:
/// SAP(i) = Σ{j ∈ neighbors(i, R)} [ Hydrophobicity(j) × (SASA(j) / SASA_max(j)) ]
///
/// Where:
/// - Neighbors are atoms/residues within radius R of atom/residue i
/// - Hydrophobicity uses the Black & Mould (1991) scale, normalized so glycine = 0
/// - SASA is the side-chain solvent accessible surface area
/// - SASA_max is the maximum SASA for that residue type
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     level (str, optional): Aggregation level for SAP calculation. Options:
///         - "atom": Calculate SAP for each atom
///         - "residue": Aggregate SAP by residue (default)
///     probe_radius (float, optional): Probe radius in Ångströms for SASA calculation. Defaults to 1.4.
///     n_points (int, optional): Number of points for SASA surface calculation. Defaults to 100.
///     model_num (int, optional): Model number to analyze (0 for first model). Defaults to 0.
///     sap_radius (float, optional): Radius in Ångströms for neighbor search. Defaults to 5.0.
///     num_threads (int, optional): Number of threads for parallel processing (0 for all cores). Defaults to 1.
///     chains (str, optional): Comma-separated chain IDs to include (e.g., "H,L").
///         If empty, includes all chains. Defaults to "".
///
/// Returns:
///     polars.DataFrame: A DataFrame with SAP scores. Columns depend on the level:
///         - atom: chain, resn, resi, insertion, atomn, atomi, sasa, sap_score
///         - residue: chain, resn, resi, insertion, sasa, sap_score
///
/// Example:
///     >>> import arpeggia
///     >>> # Residue-level SAP scores for all chains
///     >>> residue_sap = arpeggia.sap_score("structure.pdb")
///     >>> print(f"Calculated SAP for {len(residue_sap)} residues")
///     >>>
///     >>> # SAP scores for only antibody heavy and light chains
///     >>> sap_hl = arpeggia.sap_score("antibody.pdb", chains="H,L")
///     >>> print(f"Calculated SAP for {len(sap_hl)} residues")
#[pyfunction]
#[pyo3(signature = (input_file, level="residue", probe_radius=1.4, n_points=100, model_num=0, sap_radius=5.0, num_threads=1, chains=""))]
fn sap_score(
    input_file: String,
    level: &str,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
    sap_radius: f32,
    num_threads: usize,
    chains: &str,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Convert num_threads for rust-sasa
    let threads = get_num_threads(num_threads);

    // Get SAP scores based on level
    let df = match level.to_lowercase().as_str() {
        "atom" => crate::get_per_atom_sap_score(&pdb, probe_radius, n_points, model_num, sap_radius, threads, chains),
        "residue" => crate::get_per_residue_sap_score(&pdb, probe_radius, n_points, model_num, sap_radius, threads, chains),
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                "Invalid level '{}'. Must be one of: 'atom', 'residue'",
                level
            )));
        }
    };

    // Convert to PyDataFrame for Python
    Ok(PyDataFrame(df))
}

/// Calculate Shape Complementarity (SC) between two chain groups.
///
/// Shape complementarity measures how well two molecular surfaces fit together,
/// following Lawrence & Colman (1993) "Shape Complementarity at Protein/Protein Interfaces".
/// Higher SC values (closer to 1.0) indicate better geometric fit between surfaces.
/// Typical protein-protein interfaces have SC values between 0.5 and 0.7.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     groups (str): Chain groups specification, e.g., "H,L/A" for chains H,L vs chain A.
///         Both groups must be specified separated by "/".
///     model_num (int, optional): Model number to analyze (0 for first model). Defaults to 0.
///
/// Returns:
///     float: The shape complementarity score (0-1), or -1.0 if calculation fails.
///
/// Example:
///     >>> import arpeggia
///     >>> sc = arpeggia.sc("antibody_antigen.pdb", groups="H,L/A")
///     >>> print(f"SC Score: {sc:.3f}")
#[pyfunction]
#[pyo3(signature = (input_file, groups, model_num=0))]
fn sc(
    input_file: String,
    groups: &str,
    #[allow(unused_variables)]
    model_num: usize, // Kept for API compatibility; SC calculation uses all models
) -> PyResult<f64> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Calculate SC
    crate::get_sc(&pdb, groups)
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("SC calculation failed: {}", e)))
}

/// Python module for protein structure analysis.
///
/// This module provides functions for analyzing protein structures from PDB and mmCIF files,
/// including contact detection, SASA calculation, SAP score calculation, and sequence extraction.
#[pymodule]
fn arpeggia(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(contacts, m)?)?;
    m.add_function(wrap_pyfunction!(sasa, m)?)?;
    m.add_function(wrap_pyfunction!(dsasa, m)?)?;
    m.add_function(wrap_pyfunction!(relative_sasa, m)?)?;
    m.add_function(wrap_pyfunction!(sap_score, m)?)?;
    m.add_function(wrap_pyfunction!(sc, m)?)?;
    m.add_function(wrap_pyfunction!(pdb2seq, m)?)?;
    Ok(())
}
