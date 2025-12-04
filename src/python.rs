//! Python bindings for the arpeggia library using PyO3.
//!
//! This module provides Python-friendly wrappers around the core Rust functions,
//! converting Polars DataFrames to Python using pyo3-polars for efficient zero-copy data transfer.

use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

/// Load a PDB or mmCIF file and calculate atomic and ring contacts.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     groups (str, optional): Chain groups specification. Defaults to "/" (all-to-all).
///         Examples: "A,B/C,D" for chains A,B vs C,D; "A/" for chain A vs all others.
///     vdw_comp (float, optional): VdW radii compensation factor. Defaults to 0.1.
///     dist_cutoff (float, optional): Distance cutoff for neighbor searches in Ångströms. Defaults to 6.5.
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
#[pyo3(signature = (input_file, groups="/", vdw_comp=0.1, dist_cutoff=6.5))]
fn contacts(
    input_file: String,
    groups: &str,
    vdw_comp: f64,
    dist_cutoff: f64,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Get contacts
    let df = crate::get_contacts(&pdb, groups, vdw_comp, dist_cutoff);

    // Convert to PyDataFrame for Python
    Ok(PyDataFrame(df))
}

/// Load a PDB or mmCIF file and calculate solvent accessible surface area (SASA) for each atom.
///
/// Args:
///     input_file (str): Path to the PDB or mmCIF file
///     probe_radius (float, optional): Probe radius in Ångströms. Defaults to 1.4.
///     n_points (int, optional): Number of points for surface calculation. Defaults to 100.
///     model_num (int, optional): Model number to analyze (0 for first model). Defaults to 0.
///
/// Returns:
///     polars.DataFrame: A DataFrame with SASA values for each atom with columns:
///         - atomi, sasa
///         - chain, resn, resi, insertion, altloc, atomn
///
/// Example:
///     >>> import arpeggia
///     >>> sasa = arpeggia.sasa("structure.pdb", probe_radius=1.4, n_points=100)
///     >>> print(f"Calculated SASA for {len(sasa)} atoms")
#[pyfunction]
#[pyo3(signature = (input_file, probe_radius=1.4, n_points=100, model_num=0))]
fn sasa(
    input_file: String,
    probe_radius: f32,
    n_points: usize,
    model_num: usize,
) -> PyResult<PyDataFrame> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Get SASA
    let df = crate::get_atom_sasa(&pdb, probe_radius, n_points, model_num);

    // Convert to PyDataFrame for Python
    Ok(PyDataFrame(df))
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
fn sequences(input_file: String) -> PyResult<std::collections::HashMap<String, String>> {
    // Load the PDB file
    let (pdb, _warnings) = crate::load_model(&input_file);

    // Get sequences
    let seqs = crate::get_sequences(&pdb);

    Ok(seqs)
}

/// Python module for protein structure analysis.
///
/// This module provides functions for analyzing protein structures from PDB and mmCIF files,
/// including contact detection, SASA calculation, and sequence extraction.
#[pymodule]
fn arpeggia(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(contacts, m)?)?;
    m.add_function(wrap_pyfunction!(sasa, m)?)?;
    m.add_function(wrap_pyfunction!(sequences, m)?)?;
    Ok(())
}
