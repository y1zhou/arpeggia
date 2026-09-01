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
