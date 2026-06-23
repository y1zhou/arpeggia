//! Shared structure loading and preparation for calculations.

use crate::contacts::residues::ResidueExt;
use pdbtbx::*;
use std::collections::HashSet;

/// Common residue names for solvent molecules.
const SOLVENT_RESIDUES: &[&str] = &["HOH", "H2O", "D2O", "WAT", "TIP", "TIP3", "TIP4", "SPC"];

/// Common residue names for ions and terminal caps.
const ION_RESIDUES: &[&str] = &[
    "NA", "CL", "K", "CA", "MG", "ZN", "FE", "MN", "CU", "CO", "NI", "CD", "SO4", "PO4", "NO3",
    "ACE", "NH2",
];

/// Open an atomic data file with [`pdbtbx::open`] and remove non-protein residues.
pub fn load_model(input_file: &str) -> (PDB, Vec<PDBError>) {
    let (mut pdb, errors) = pdbtbx::ReadOptions::default()
        .set_only_atomic_coords(true)
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(input_file)
        .unwrap();

    pdb.remove_residues_by(|res| res.resn().is_none());

    (pdb, errors)
}

/// Parse the chain groups from the input string.
///
/// Only checks the first two fields separated by `/`. If one of the groups is
/// unspecified, all remaining chains from `all_chains` are used.
///
/// # Panics
/// This function will panic if the input format is invalid or if one of the
/// groups ends up empty.
pub fn parse_groups(
    all_chains: &HashSet<String>,
    groups: &str,
) -> (HashSet<String>, HashSet<String>) {
    let sel_vec: Vec<&str> = groups.split('/').collect();
    assert!(
        (sel_vec.len() >= 2),
        "Invalid chain groups format! Use '/' for all-to-all comparisons."
    );
    let ligand_chains = sel_vec.first().unwrap_or(&"");
    let receptor_chains = sel_vec.get(1).unwrap_or(&"");

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

    if ligand.is_empty() && receptor.is_empty() {
        return (all_chains.clone(), all_chains.clone());
    }

    if ligand.is_empty() {
        ligand = all_chains.difference(&receptor).cloned().collect();
    } else if receptor.is_empty() {
        receptor = all_chains.difference(&ligand).cloned().collect();
    }

    assert!(
        !(ligand.is_empty() || receptor.is_empty()),
        "Empty chain groups!"
    );

    (ligand, receptor)
}

/// Structure preparation policy shared by SASA, SAP, dSASA, and SC.
#[derive(Clone, Debug, Default)]
pub(crate) struct StructurePreparation {
    model_num: usize,
    remove_hydrogens: bool,
    remove_solvent_and_ions: bool,
    chains: HashSet<String>,
}

impl StructurePreparation {
    pub(crate) fn new(model_num: usize) -> Self {
        Self {
            model_num,
            ..Self::default()
        }
    }

    pub(crate) fn remove_hydrogens(mut self, remove_hydrogens: bool) -> Self {
        self.remove_hydrogens = remove_hydrogens;
        self
    }

    pub(crate) fn remove_solvent_and_ions(mut self, remove_solvent_and_ions: bool) -> Self {
        self.remove_solvent_and_ions = remove_solvent_and_ions;
        self
    }

    pub(crate) fn chain_csv(mut self, chains: &str) -> Self {
        self.chains = parse_chain_string(chains);
        self
    }

    pub(crate) fn chain_set(mut self, chains: &HashSet<String>) -> Self {
        self.chains = chains.clone();
        self
    }

    pub(crate) fn prepare(&self, pdb: &PDB) -> PDB {
        let mut pdb_prepared = filter_pdb_by_model(pdb, self.model_num);
        assert!(
            pdb_prepared.model_count() > 0,
            "No models exist after preparation step; check `model_num`"
        );

        if !self.chains.is_empty() {
            pdb_prepared.remove_chains_by(|chain| !self.chains.contains(chain.id()));
            assert!(
                pdb_prepared.chain_count() > 0,
                "No chains exist after preparation step; check `chains`"
            );
        }

        if self.remove_hydrogens {
            pdb_prepared.remove_atoms_by(|atom| atom.element() == Some(&Element::H));
        }

        if self.remove_solvent_and_ions {
            pdb_prepared.remove_residues_by(|residue| {
                let resn = residue.name().unwrap_or("");
                SOLVENT_RESIDUES.contains(&resn) || ION_RESIDUES.contains(&resn)
            });
        }

        pdb_prepared
    }
}

/// Prepare a structure for surface-based calculations.
pub(crate) fn prepare_structure(
    pdb: &PDB,
    model_num: usize,
    remove_hydrogens: bool,
    chains: &str,
) -> PDB {
    StructurePreparation::new(model_num)
        .remove_hydrogens(remove_hydrogens)
        .remove_solvent_and_ions(true)
        .chain_csv(chains)
        .prepare(pdb)
}

/// Prepare a structure for a specific set of chain IDs.
pub(crate) fn prepare_structure_with_chains(
    pdb: &PDB,
    model_num: usize,
    remove_hydrogens: bool,
    chains: &HashSet<String>,
) -> PDB {
    StructurePreparation::new(model_num)
        .remove_hydrogens(remove_hydrogens)
        .remove_solvent_and_ions(true)
        .chain_set(chains)
        .prepare(pdb)
}

/// Filter a PDB structure to keep only the specified model.
pub(crate) fn filter_pdb_by_model(pdb: &PDB, model_num: usize) -> PDB {
    let all_model_nums: Vec<usize> = pdb.models().map(|m| m.serial_number()).collect();
    let picked_model_num = if model_num == 0 {
        *all_model_nums
            .first()
            .expect("No models found in the structure")
    } else {
        if !all_model_nums.contains(&model_num) {
            panic!("Model {model_num} not found in the structure");
        }
        model_num
    };

    let picked_model_idx = all_model_nums
        .into_iter()
        .position(|m| m == picked_model_num)
        .unwrap_or(0);

    let mut pdb_filtered = pdb.clone();
    pdb_filtered.remove_models_except(&[picked_model_idx]);
    pdb_filtered
}

fn parse_chain_string(chains: &str) -> HashSet<String> {
    if chains.is_empty() {
        HashSet::new()
    } else {
        chains
            .split(',')
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn good_group_splits() {
        let chains: HashSet<String> = HashSet::from(["A", "B", "C", "D"].map(|c| c.to_string()));

        assert_eq!(
            (
                HashSet::from(["A", "B"].map(|c| c.to_string())),
                HashSet::from(["C", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "A,B/C,D")
        );

        assert_eq!(
            (
                HashSet::from(["A"].map(|c| c.to_string())),
                HashSet::from(["C", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "A/C,D")
        );

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
        let chains: HashSet<String> = HashSet::from(["A", "B", "C", "D"].map(|c| c.to_string()));
        parse_groups(&chains, "");
    }

    #[test]
    #[should_panic(expected = "Empty chain groups!")]
    fn missing_groups_in_split() {
        let chains: HashSet<String> = HashSet::from(["A", "B", "C"].map(|c| c.to_string()));
        parse_groups(&chains, "A,B,C/");
    }

    #[test]
    fn removes_zero_occupancy_only_when_requested_by_callers() {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/1ubq.pdb");

        let (mut pdb, _) = load_model(&path);
        let initial_atom_count = pdb.atom_count();

        pdb.remove_atoms_by(|atom| atom.occupancy() == 0.0);
        let final_atom_count = pdb.atom_count();

        assert_eq!(initial_atom_count, final_atom_count);
    }
}
