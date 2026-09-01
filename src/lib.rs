#![warn(missing_docs)]
#![warn(unreachable_pub)]
#![doc = include_str!("../README.md")]

//! # Arpeggia Library
//!
//! This library provides functionality for analyzing protein-protein interactions
//! in PDB and mmCIF files. It can identify various types of interactions including
//! hydrogen bonds, ionic interactions, aromatic interactions, and more.
//!
//! The library returns results as Polars `DataFrames`, which can be easily converted
//! to various output formats or used directly in Python via `PyO3` bindings.

mod clustering;
mod contacts;
mod diagnostics;
mod metadata;
mod pairwise_rmsd;
mod rmsd;
mod sap;
mod sasa;
mod sc;
mod sequences;
mod structure;
mod utils;

// Re-export key public types
pub use clustering::{ClusterOptions, ClusteringMethod, cluster_pairwise_rmsd};
pub use contacts::ProtonationMode;
pub use diagnostics::{Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
pub use metadata::{BondEndpoint, StructureMetadata, get_seqres, read_metadata};
pub use pairwise_rmsd::{
    PairwiseRmsdMatrix, PairwiseRmsdOptions, StructureObservation, get_pairwise_rmsd,
    get_pairwise_rmsd_matrix, read_pairwise_matrix, read_structure_observations,
};
pub use rmsd::{AtomSubset, get_rmsd, kabsch_rmsd, validate_rmsd_selections};
pub use structure::load_model;
pub use utils::run_with_threads;

// Re-export public functions from modules
pub use contacts::analyze_contacts;
pub use sap::{get_per_atom_sap_score, get_per_residue_sap_score};
pub use sasa::{
    DsasaResult, get_atom_sasa, get_chain_sasa, get_dsasa_components, get_max_asa,
    get_relative_sasa, get_residue_sasa,
};
pub use sc::{ScResult, get_sc_details};
pub use sequences::get_sequences;

// Python bindings module (only compiled when python feature is enabled)
#[cfg(feature = "python")]
mod python;

// #[cfg(feature = "python")]
// pub use python::*;
