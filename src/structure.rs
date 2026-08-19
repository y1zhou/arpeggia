//! Shared structure loading and preparation for calculations.

use crate::contacts::residues::ResidueExt;
use crate::{Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
use pdbtbx::*;
use std::collections::HashSet;

/// Common residue names for solvent molecules.
const SOLVENT_RESIDUES: &[&str] = &["HOH", "H2O", "D2O", "WAT", "TIP", "TIP3", "TIP4", "SPC"];

/// Common residue names for ions and terminal caps.
const ION_RESIDUES: &[&str] = &[
    "NA", "CL", "K", "CA", "MG", "ZN", "FE", "MN", "CU", "CO", "NI", "CD", "SO4", "PO4", "NO3",
    "ACE", "NH2",
];

/// Open an atomic data file and remove unsupported residues.
pub fn load_model(input_file: &str) -> ArpeggiaResult<Analysis<PDB>> {
    if !std::path::Path::new(input_file).try_exists()? {
        return Err(ArpeggiaError::Io(std::io::Error::new(
            std::io::ErrorKind::NotFound,
            format!("input file does not exist: {input_file}"),
        )));
    }
    let (mut pdb, diagnostics) = pdbtbx::ReadOptions::default()
        .set_only_atomic_coords(true)
        .set_level(pdbtbx::StrictnessLevel::Loose)
        .read(input_file)
        .map_err(|errors| ArpeggiaError::Parse(format_pdb_errors(&errors)))?;

    let invalidating: Vec<_> = diagnostics
        .iter()
        .filter(|error| {
            matches!(
                error.level(),
                ErrorLevel::BreakingError | ErrorLevel::InvalidatingError
            )
        })
        .collect();
    if !invalidating.is_empty() {
        return Err(ArpeggiaError::Parse(
            invalidating
                .into_iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join("; "),
        ));
    }

    pdb.remove_residues_by(|res| res.resn().is_none());

    let mut warnings = diagnostics
        .into_iter()
        .map(|error| AnalysisWarning::new(WarningCode::Parser, error.to_string()))
        .collect::<Vec<_>>();
    warnings.extend(select_conformers(&mut pdb));
    Ok(Analysis::new(pdb, warnings))
}

fn format_pdb_errors(errors: &[PDBError]) -> String {
    errors
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join("; ")
}

pub(crate) fn select_conformers(pdb: &mut PDB) -> Vec<AnalysisWarning> {
    let mut warnings = Vec::new();
    for model in pdb.models_mut() {
        let model_number = model.serial_number();
        for chain in model.chains_mut() {
            let chain_id = chain.id().to_string();
            for residue in chain.residues_mut() {
                let choices = residue
                    .conformers()
                    .filter_map(|conformer| {
                        let location = conformer.alternative_location()?;
                        let atom_count = conformer.atom_count();
                        (atom_count > 0).then(|| {
                            let mean_occupancy =
                                conformer.atoms().map(Atom::occupancy).sum::<f64>()
                                    / atom_count as f64;
                            (location.to_string(), mean_occupancy)
                        })
                    })
                    .collect::<Vec<_>>();
                if choices.len() < 2 {
                    continue;
                }
                // [WARNING] Selecting one highest-occupancy conformer is deterministic but does
                // not model conformational ensembles or their population-weighted observables.
                let Some(selected) = choices
                    .iter()
                    .max_by(
                        |(left_name, left_occupancy), (right_name, right_occupancy)| {
                            left_occupancy.total_cmp(right_occupancy).then_with(|| {
                                conformer_tie_rank(right_name).cmp(&conformer_tie_rank(left_name))
                            })
                        },
                    )
                    .map(|(name, _)| name.clone())
                else {
                    continue;
                };
                residue.remove_conformers_by(|conformer| {
                    conformer
                        .alternative_location()
                        .is_some_and(|location| location != selected)
                });
                warnings.push(AnalysisWarning::new(
                    WarningCode::ConformerSelected,
                    format!(
                        "selected alternate conformer {selected} for model {model_number}, chain \
                         {chain_id}, residue {}{}",
                        residue.serial_number(),
                        residue.insertion_code().unwrap_or("")
                    ),
                ));
            }
        }
    }
    warnings
}

fn conformer_tie_rank(name: &str) -> (bool, &str) {
    (name != "A", name)
}

/// Parse the chain groups from the input string.
///
/// Exactly one `/` is required. If one side is unspecified, all remaining
/// chains from `all_chains` are used. Invalid or empty selections return a
/// typed argument error.
pub fn parse_groups(
    all_chains: &HashSet<String>,
    groups: &str,
) -> ArpeggiaResult<(HashSet<String>, HashSet<String>)> {
    let sel_vec: Vec<&str> = groups.split('/').collect();
    if sel_vec.len() != 2 {
        return Err(ArpeggiaError::InvalidArgument(
            "chain groups must contain exactly one '/'; use '/' for all-to-all comparisons".into(),
        ));
    }
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
        return Ok((all_chains.clone(), all_chains.clone()));
    }

    if ligand.is_empty() {
        ligand = all_chains.difference(&receptor).cloned().collect();
    } else if receptor.is_empty() {
        receptor = all_chains.difference(&ligand).cloned().collect();
    }

    if ligand.is_empty() || receptor.is_empty() {
        return Err(ArpeggiaError::InvalidArgument(
            "chain selection produced an empty group".into(),
        ));
    }
    let unknown = ligand
        .union(&receptor)
        .filter(|chain| !all_chains.contains(*chain))
        .cloned()
        .collect::<Vec<_>>();
    if !unknown.is_empty() {
        return Err(ArpeggiaError::InvalidArgument(format!(
            "unknown chain identifiers: {}",
            unknown.join(",")
        )));
    }

    Ok((ligand, receptor))
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
        select_conformers(&mut pdb_prepared);
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

pub(crate) fn selected_model(pdb: &PDB, model_num: usize) -> ArpeggiaResult<&Model> {
    if model_num == 0 {
        pdb.models().next()
    } else {
        pdb.models()
            .find(|model| model.serial_number() == model_num)
    }
    .ok_or_else(|| ArpeggiaError::InvalidArgument(format!("model {model_num} does not exist")))
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
    fn missing_input_is_a_typed_error() {
        let error = load_model("this-file-does-not-exist.pdb").unwrap_err();
        assert!(matches!(error, ArpeggiaError::Io(_)));
    }

    #[test]
    fn load_selects_highest_occupancy_conformer_with_a_tie_break() {
        let input =
            b"ATOM      1  CB AALA A   1       0.000   0.000   0.000  0.50 20.00           C  \n\
ATOM      2  CB BALA A   1       1.000   0.000   0.000  0.50 20.00           C  \n\
ATOM      3  CB AALA A   2       2.000   0.000   0.000  0.30 20.00           C  \n\
ATOM      4  CB BALA A   2       3.000   0.000   0.000  0.70 20.00           C  \n\
END                                                                             \n";
        let path =
            std::env::temp_dir().join(format!("arpeggia-conformers-{}.pdb", std::process::id()));
        std::fs::write(&path, input).unwrap();
        let analysis = load_model(path.to_str().unwrap()).unwrap();
        std::fs::remove_file(path).unwrap();

        let altlocs = analysis
            .value
            .atoms_with_hierarchy()
            .map(|entity| {
                entity
                    .conformer()
                    .alternative_location()
                    .unwrap()
                    .to_string()
            })
            .collect::<Vec<_>>();
        assert_eq!(altlocs, ["A", "B"]);
        assert_eq!(
            analysis
                .warnings
                .iter()
                .filter(|warning| warning.code == WarningCode::ConformerSelected)
                .count(),
            2
        );
    }

    #[test]
    fn good_group_splits() {
        let chains: HashSet<String> = HashSet::from(["A", "B", "C", "D"].map(|c| c.to_string()));

        assert_eq!(
            (
                HashSet::from(["A", "B"].map(|c| c.to_string())),
                HashSet::from(["C", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "A,B/C,D").unwrap()
        );

        assert_eq!(
            (
                HashSet::from(["A"].map(|c| c.to_string())),
                HashSet::from(["C", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "A/C,D").unwrap()
        );

        assert_eq!(
            (
                HashSet::from(["A", "B"].map(|c| c.to_string())),
                HashSet::from(["C", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "/C,D").unwrap()
        );
        assert_eq!(
            (
                HashSet::from(["C"].map(|c| c.to_string())),
                HashSet::from(["A", "B", "D"].map(|c| c.to_string()))
            ),
            parse_groups(&chains, "C/").unwrap()
        );
        assert_eq!(
            (chains.clone(), chains.clone()),
            parse_groups(&chains, "/").unwrap()
        );
    }

    #[test]
    fn empty_group_splits() {
        let chains: HashSet<String> = HashSet::from(["A", "B", "C", "D"].map(|c| c.to_string()));
        assert!(parse_groups(&chains, "").is_err());
    }

    #[test]
    fn missing_groups_in_split() {
        let chains: HashSet<String> = HashSet::from(["A", "B", "C"].map(|c| c.to_string()));
        assert!(parse_groups(&chains, "A,B,C/").is_err());
    }

    #[test]
    fn removes_zero_occupancy_only_when_requested_by_callers() {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/1ubq.pdb");

        let mut pdb = load_model(&path).unwrap().value;
        let initial_atom_count = pdb.atom_count();

        pdb.remove_atoms_by(|atom| atom.occupancy() == 0.0);
        let final_atom_count = pdb.atom_count();

        assert_eq!(initial_atom_count, final_atom_count);
    }
}
