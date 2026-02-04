#![warn(missing_docs)]
#![doc = include_str!("../README.md")]

//! # Arpeggia Library
//!
//! This library provides functionality for analyzing protein-protein interactions
//! in PDB and mmCIF files. It can identify various types of interactions including
//! hydrogen bonds, ionic interactions, aromatic interactions, and more.
//!
//! The library returns results as Polars `DataFrames`, which can be easily converted
//! to various output formats or used directly in Python via `PyO3` bindings.

mod chains;
mod contacts;
mod interactions;
mod residues;
mod sap;
mod sasa;
mod sc;
mod sequences;
mod utils;

// Re-export key public types
pub use interactions::{InteractingEntity, Interaction, Interactions, ResultEntry};
pub use residues::{Plane, ResidueExt, ResidueId};
pub use utils::{
    DataFrameFileType, get_num_threads, load_model, parse_groups, sum_float_col, write_df_to_file,
};

// Re-export public functions from modules
pub use contacts::get_contacts;
pub use sap::{get_per_atom_sap_score, get_per_residue_sap_score};
pub use sasa::{
    get_atom_sasa, get_chain_sasa, get_dsasa, get_max_asa, get_relative_sasa, get_residue_sasa,
};
pub use sc::get_sc;
pub use sequences::get_sequences;

// Python bindings module (only compiled when python feature is enabled)
#[cfg(feature = "python")]
mod python;

// #[cfg(feature = "python")]
// pub use python::*;
