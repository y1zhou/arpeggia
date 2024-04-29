use std::collections::HashSet;

use crate::residues::ResidueExt;
use pdbtbx::*;

pub fn load_model(input_file: &String) -> (PDB, Vec<PDBError>) {
    // Load file as complex structure
    let (mut pdb, errors) = pdbtbx::open(input_file, StrictnessLevel::Loose).unwrap();

    // Remove non-protein residues from model
    pdb.remove_residues_by(|res| res.resn().is_none());

    (pdb, errors)
}

pub fn parse_groups(
    all_chains: &HashSet<String>,
    groups: &str,
) -> (HashSet<String>, HashSet<String>) {
    // Parse the first two fields in groups
    let sel_vec: Vec<&str> = groups.split('/').collect();
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
