//! Sequence extraction from protein structures.
//!
//! This module provides functions for extracting amino acid sequences
//! from PDB structures.

use crate::contacts::chains::ChainExt;
use pdbtbx::*;
use std::collections::HashMap;

/// Get sequences of all chains in a PDB structure.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
///
/// # Returns
///
/// A `HashMap` mapping chain IDs to their sequences as strings.
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
    pdb.chains()
        .map(|chain| (chain.id().to_string(), chain.pdb_seq().join("")))
        .collect()
}
