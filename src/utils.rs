use crate::residues::ResidueExt;
use pdbtbx::*;
use polars::prelude::*;
use std::{collections::HashSet, path::Path};

/// Convert thread count from usize to isize for rust-sasa.
///
/// The rust-sasa library uses isize for thread count where:
/// - -1 means use all available cores
/// - 1 means single-threaded
/// - n > 1 means use n threads
///
/// This function converts from the user-facing usize convention where:
/// - 0 means use all available cores
/// - n > 0 means use n threads
///
/// # Arguments
///
/// * `num_threads` - Number of threads (0 for all cores)
///
/// # Returns
///
/// Thread count as isize suitable for rust-sasa
pub fn threads_to_isize(num_threads: usize) -> isize {
    if num_threads == 0 {
        -1
    } else {
        num_threads as isize
    }
}

/// Sum the SASA column of a DataFrame using polars lazy aggregation.
///
/// # Arguments
///
/// * `df` - DataFrame containing a "sasa" column
///
/// # Returns
///
/// The sum of all SASA values, or 0.0 if the column is empty or doesn't exist.
pub fn sum_sasa(df: &DataFrame) -> f32 {
    df.clone()
        .lazy()
        .select([col("sasa").sum()])
        .collect()
        .unwrap()
        .column("sasa")
        .unwrap()
        .f32()
        .unwrap()
        .get(0)
        .unwrap_or(0.0)
}

/// Open an atomic data file with [`pdbtbx::open`] and remove non-protein residues.
pub fn load_model(input_file: &String) -> (PDB, Vec<PDBError>) {
    // Load file as complex structure
    let (mut pdb, errors) = pdbtbx::ReadOptions::default()
        .set_only_atomic_coords(true)
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(input_file)
        .unwrap();

    // Remove non-protein residues from model
    pdb.remove_residues_by(|res| res.resn().is_none());

    (pdb, errors)
}

/// Parse the chain groups from the input string.
/// Only checks the first two fields separated by `/`.
/// If one of the groups is unspecified, all remaining chains from `all_chains` are used.
pub fn parse_groups(
    all_chains: &HashSet<String>,
    groups: &str,
) -> (HashSet<String>, HashSet<String>) {
    // Parse the first two fields in groups
    let sel_vec: Vec<&str> = groups.split('/').collect();
    if sel_vec.len() < 2 {
        panic!("Invalid chain groups format! Use '/' for all-to-all comparisons.")
    }
    let ligand_chains = sel_vec.first().unwrap_or(&"");
    let receptor_chains = sel_vec.get(1).unwrap_or(&"");

    // Create a HashSet of chains for the ligand and receptor
    let mut ligand: HashSet<String> = ligand_chains
        .split(',')
        .map(|c| c.to_string())
        .filter(|c| !c.is_empty())
        .collect();
    let mut receptor: HashSet<String> = receptor_chains
        .split(',')
        .map(|c| c.to_string())
        .filter(|c| !c.is_empty())
        .collect();

    // If both groups are empty, perform all-to-all comparisons
    if ligand.is_empty() && receptor.is_empty() {
        return (all_chains.clone(), all_chains.clone());
    }

    // If there are no ligand or receptor chains, use all remaining chains
    if ligand.is_empty() {
        ligand = all_chains.difference(&receptor).cloned().collect();
    } else if receptor.is_empty() {
        receptor = all_chains.difference(&ligand).cloned().collect();
    }

    // Panic if there are no chains in one of the groups
    if ligand.is_empty() || receptor.is_empty() {
        panic!("Empty chain groups!")
    }

    (ligand, receptor)
}

/// Write a DataFrame to a CSV file
pub fn write_df_to_file(df: &mut DataFrame, file_path: &Path, file_type: DataFrameFileType) {
    let file_suffix = file_type.to_string();
    let mut file = std::fs::File::create(file_path.with_extension(file_suffix)).unwrap();
    match file_type {
        DataFrameFileType::Csv => {
            CsvWriter::new(&mut file).finish(df).unwrap();
        }
        DataFrameFileType::Parquet => {
            ParquetWriter::new(&mut file).finish(df).unwrap();
        }
        DataFrameFileType::Json => {
            JsonWriter::new(&mut file)
                .with_json_format(JsonFormat::Json)
                .finish(df)
                .unwrap();
        }
        DataFrameFileType::NDJson => {
            JsonWriter::new(&mut file)
                .with_json_format(JsonFormat::JsonLines)
                .finish(df)
                .unwrap();
        }
    }
}

/// File format for writing DataFrames.
#[derive(clap::ValueEnum, Clone, Debug, Copy)]
pub enum DataFrameFileType {
    /// Comma-separated values
    Csv,
    /// Parquet columnar storage
    Parquet,
    /// Standard JSON
    Json,
    /// Newline-delimited JSON
    NDJson,
}

impl std::fmt::Display for DataFrameFileType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            DataFrameFileType::Csv => write!(f, "csv"),
            DataFrameFileType::Parquet => write!(f, "parquet"),
            DataFrameFileType::Json => write!(f, "json"),
            DataFrameFileType::NDJson => write!(f, "ndjson"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn good_group_splits() {
        let chains: HashSet<String> = HashSet::from(["A", "B", "C", "D"].map(|c| c.to_string()));

        // All chains are specified
        assert_eq!(
            (
                HashSet::from(["A", "B"].map(|c| c.to_string())),
                HashSet::from(["C", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "A,B/C,D")
        );

        // One chain missing (ignored in searches)
        assert_eq!(
            (
                HashSet::from(["A"].map(|c| c.to_string())),
                HashSet::from(["C", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "A/C,D")
        );

        // Empty group
        assert_eq!(
            (
                HashSet::from(["A", "B"].map(|c| c.to_string())),
                HashSet::from(["C", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "/C,D")
        );
        assert_eq!(
            (
                HashSet::from(["C"].map(|c| c.to_string())),
                HashSet::from(["A", "B", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "C/")
        );
        assert_eq!((chains.clone(), chains.clone()), parse_groups(&chains, "/"));
    }

    #[test]
    #[should_panic(expected = "Invalid chain groups format! Use '/' for all-to-all comparisons.")]
    fn empty_group_splits() {
        // Nothing passed
        let chains: HashSet<String> = HashSet::from(["A", "B", "C", "D"].map(|c| c.to_string()));
        parse_groups(&chains, "");
    }

    #[test]
    #[should_panic(expected = "Empty chain groups!")]
    fn missing_groups_in_split() {
        // Nothing passed
        let chains: HashSet<String> = HashSet::from(["A", "B", "C"].map(|c| c.to_string()));
        parse_groups(&chains, "A,B,C/");
    }

    #[test]
    fn test_remove_atoms_zero_occupancy() {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/1ubq.pdb");

        let (mut pdb, _) = load_model(&path);
        let initial_atom_count = pdb.atom_count();

        // Remove atoms with zero occupancy (1ubq.pdb has all atoms with occupancy 1.0)
        pdb.remove_atoms_by(|atom| atom.occupancy() == 0.0);
        let final_atom_count = pdb.atom_count();

        // No atoms should be removed since all have occupancy 1.0
        assert_eq!(
            initial_atom_count, final_atom_count,
            "Atom count should remain the same since all atoms have occupancy 1.0"
        );
    }
}
