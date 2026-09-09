pub(crate) mod align_seqs;
pub(crate) mod cluster_structs;
pub(crate) mod contacts;
pub(crate) mod dsasa;
mod output;
pub(crate) mod pdb2seq;
pub(crate) mod relative_sasa;
pub(crate) mod rmsd;
pub(crate) mod sap;
pub(crate) mod sasa;
pub(crate) mod sc;

use arpeggia::{ArpeggiaError, ArpeggiaResult, load_model};
use output::{DataFrameFileType, prepare_df_output_dir, write_df_to_file, write_df_to_new_file};
use pdbtbx::PDB;
use std::path::Path;
use tracing::warn;

fn load_input(path: &Path) -> ArpeggiaResult<PDB> {
    let path = path.canonicalize().map_err(|error| {
        ArpeggiaError::Io(std::io::Error::new(
            error.kind(),
            format!("cannot resolve input {}: {error}", path.display()),
        ))
    })?;
    let input = path
        .to_str()
        .ok_or_else(|| ArpeggiaError::InvalidArgument("input path is not valid UTF-8".into()))?;
    let analysis = load_model(input)?;
    for warning in analysis.warnings {
        warn!("{warning}");
    }
    Ok(analysis.value)
}

// Shared by RMSD and clustering so selection syntax is available in --help.
pub(crate) const RESIDUE_SELECTION_HELP: &str = "Residue selections:
  Empty selects all eligible residues, independently for fitting and evaluation.
  A                 all residues in chain A
  A:1-100,A:110-120,B inclusive author-number ranges plus all of chain B
  A:-5--1,B:10A-20   negative numbers and insertion codes
  A:10              all insertion variants at author number 10
  A:10A             only insertion 10A
Repeat the chain in every comma-separated clause (A:1,A:3).
An empty --rmsd-residues does not inherit --superpose-residues.
Selections apply to both structures; rmsd --align-seqs uses reference numbering.
Fitting needs at least three non-collinear atom pairs; evaluation needs one.
Examples, schemas, and cache behavior:
https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md";
