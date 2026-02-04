use crate::residues::ResidueExt;
use pdbtbx::*;
use polars::prelude::*;
use std::{collections::HashSet, path::Path};

/// Execute a parallel operation with the configured thread limit.
/// Uses rayon's thread pool with the specified number of threads.
pub fn run_with_threads<F, T>(num_threads: isize, f: F) -> T
where
    F: FnOnce() -> T + Send,
    T: Send,
{
    let n_threads = num_threads.max(0) as usize;
    if n_threads == 0 || n_threads >= rayon::current_num_threads() {
        // Use global pool directly when auto or enough threads
        f()
    } else {
        // Create a scoped pool with limited threads
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build()
            .unwrap_or_else(|_| {
                // Fallback to global pool if creation fails
                rayon::ThreadPoolBuilder::new()
                    .build()
                    .expect("Failed to create thread pool")
            });
        pool.install(f)
    }
}

/// Sum the SASA column of a `DataFrame`.
///
/// # Arguments
///
/// * `df` - `DataFrame` containing a "sasa" column
///
/// # Returns
///
/// The sum of all SASA values, or 0.0 if the column is empty or doesn't exist.
pub fn sum_float_col(df: &DataFrame, colname: &str) -> f32 {
    df.column(colname)
        .unwrap()
        .f32()
        .unwrap()
        .sum()
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
///
/// # Panics
/// This function will panic if the input format is invalid or if one of the groups ends up empty.
pub fn parse_groups(
    all_chains: &HashSet<String>,
    groups: &str,
) -> (HashSet<String>, HashSet<String>) {
    // Parse the first two fields in groups
    let sel_vec: Vec<&str> = groups.split('/').collect();
    assert!(
        (sel_vec.len() >= 2),
        "Invalid chain groups format! Use '/' for all-to-all comparisons."
    );
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
    assert!(
        !(ligand.is_empty() || receptor.is_empty()),
        "Empty chain groups!"
    );

    (ligand, receptor)
}

/// Write a `DataFrame` to a CSV file
///
/// # Panics
/// This function will panic if the file cannot be created or written to.
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

/// File format for writing `DataFrames`.
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
