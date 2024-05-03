use std::collections::HashSet;

use crate::residues::ResidueExt;
use pdbtbx::*;

/// Open an atomic data file with [`pdbtbx::open`] and remove non-protein residues.
pub fn load_model(input_file: &String) -> (PDB, Vec<PDBError>) {
    // Load file as complex structure
    let (mut pdb, errors) = pdbtbx::open(input_file, StrictnessLevel::Loose).unwrap();

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
        panic!("Invalid chain groups format!")
    }
    let ligand_chains = sel_vec.first().unwrap();
    let receptor_chains = sel_vec.get(1).unwrap();

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

    // If there are no ligand or receptor chains, use all remaining chains
    if ligand.is_empty() ^ receptor.is_empty() {
        if ligand.is_empty() {
            ligand = all_chains.difference(&receptor).cloned().collect();
        } else {
            receptor = all_chains.difference(&ligand).cloned().collect();
        }
    }

    if ligand.is_empty() || receptor.is_empty() {
        panic!("Empty chain groups!")
    }

    (ligand, receptor)
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
    }

    #[test]
    #[should_panic(expected = "Invalid chain groups format!")]
    fn empty_group_splits() {
        // Nothing passed
        let chains: HashSet<String> = HashSet::from(["A", "B", "C", "D"].map(|c| c.to_string()));
        parse_groups(&chains, "");
    }

    #[test]
    #[should_panic(expected = "Empty chain groups!")]
    fn missing_groups_in_split() {
        // Nothing passed
        let chains: HashSet<String> = HashSet::from(["A", "B", "C", "D"].map(|c| c.to_string()));
        parse_groups(&chains, "/");
    }
}
