pub mod contacts;
pub mod dsasa;
pub mod pdb2seq;
pub mod relative_sasa;
pub mod sap;
pub mod sasa;
pub mod sc;

use arpeggia::{ArpeggiaError, ArpeggiaResult, load_model};
use pdbtbx::PDB;
use std::path::Path;
use tracing::warn;

pub(crate) fn load_input(path: &Path) -> ArpeggiaResult<PDB> {
    let path = path.canonicalize()?;
    let input = path
        .to_str()
        .ok_or_else(|| ArpeggiaError::InvalidArgument("input path is not valid UTF-8".into()))?;
    let analysis = load_model(input)?;
    for warning in analysis.warnings {
        warn!("{warning}");
    }
    Ok(analysis.value)
}
