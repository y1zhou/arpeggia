use super::{
    find_hydrogen_bond, find_hydrophobic_contact, find_ionic_bond, find_vdw_contact,
    find_weak_hydrogen_bond, Interaction, ResultEntry,
};
use crate::{
    residues::{ResidueExt, Ring},
    utils::{hierarchy_to_entity, parse_groups},
};

use pdbtbx::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// The workhorse struct for identifying interactions in the model
pub struct InteractionComplex {
    /// All information present in the atomic model
    pub model: PDB,
    /// All ligand chains
    pub ligand: HashSet<String>,
    /// All receptor chains
    pub receptor: HashSet<String>,
    /// Compensation factor for VdW radii dependent interaction types
    pub vdw_comp_factor: f64,
    /// Distance cutoff when searching for neighboring atoms
    pub interacting_threshold: f64,

    /// Maps residue names to unique indices
    res2idx: HashMap<String, HashMap<(isize, String), usize>>,
    rings: HashMap<(String, isize, String), Ring>,
}

impl InteractionComplex {
    pub fn new(model: PDB, groups: &str, vdw_comp_factor: f64, interacting_threshold: f64) -> Self {
        // Parse all chains and input chain groups
        let all_chains: HashSet<String> = model.par_chains().map(|c| c.id().to_string()).collect();
        let (ligand, receptor) = parse_groups(&all_chains, groups);

        // Build a mapping of residue names to indices
        let res2idx = build_residue_index(&model);

        // Build a mapping of ring residue names to ring centers and normals
        let rings = build_ring_positions(&model);

        Self {
            model,
            ligand,
            receptor,
            vdw_comp_factor,
            interacting_threshold,
            res2idx,
            rings,
        }
    }

    fn is_neighboring_res_pair(
        &self,
        x: &AtomConformerResidueChainModel,
        y: &AtomConformerResidueChainModel,
    ) -> bool {
        if x.chain().id() != y.chain().id() {
            return false;
        }
        let (x_resi, x_insertion) = x.residue().id();
        let x_altloc = match x_insertion {
            Some(insertion) => insertion.to_string(),
            None => "".to_string(),
        };
        let x_idx = self.res2idx[x.chain().id()][&(x_resi, x_altloc)];
        let (y_resi, y_insertion) = y.residue().id();
        let y_altloc = match y_insertion {
            Some(insertion) => insertion.to_string(),
            None => "".to_string(),
        };
        let y_idx = self.res2idx[y.chain().id()][&(y_resi, y_altloc)];

        match x_idx {
            0 => y_idx == x_idx + 1,
            _ => (y_idx == x_idx - 1) | (y_idx == x_idx + 1),
        }
    }
}

pub trait Interactions {
    /// Get all protein-protein interactions between the ligand and receptor
    fn get_ppi(&self);
    /// Get all atomic interactions between the ligand and receptor
    fn get_atomic_contacts(&self) -> Vec<ResultEntry>;
    /// Get all ring interactions between the ligand and receptor
    fn get_ring_contacts(&self);
}

impl Interactions for InteractionComplex {
    fn get_ppi(&self) {
        todo!()
    }

    fn get_atomic_contacts(&self) -> Vec<ResultEntry> {
        let tree = self.model.create_hierarchy_rtree();
        let max_radius_squared = self.interacting_threshold * self.interacting_threshold;

        // Find all atoms within the radius of the ligand atoms
        let ligand_neighbors: Vec<(
            AtomConformerResidueChainModel,
            AtomConformerResidueChainModel,
        )> = self
            .model
            .atoms_with_hierarchy()
            .filter(|x| {
                self.ligand.contains(x.chain().id()) & (x.atom().element().unwrap() != &Element::H)
            })
            .flat_map(|x| {
                let neighbor_entities: Vec<(
                    AtomConformerResidueChainModel,
                    AtomConformerResidueChainModel,
                )> = tree
                    .locate_within_distance(x.atom().pos(), max_radius_squared)
                    .filter(|y| {
                        self.receptor.contains(y.chain().id())
                            & (y.atom().element().unwrap() != &Element::H)
                            // Skip lower resi if both entities are on the same chain
                            & !((x.chain().id() == y.chain().id())
                                & (x.residue().serial_number() > y.residue().serial_number()))
                            // Skip neighboring residues
                            & !self.is_neighboring_res_pair(&x, y)
                            // Skip atoms in the same residue
                            & !((x.chain().id() == y.chain().id())
                                & (x.residue().serial_number() == y.residue().serial_number())
                                & (x.residue().insertion_code() == y.residue().insertion_code()))
                    })
                    .map(|y| (x.clone(), y.clone()))
                    .collect();

                match neighbor_entities.len() {
                    0 => vec![],
                    _ => neighbor_entities,
                }
            })
            .collect();

        ligand_neighbors
            .par_iter()
            .filter_map(|x| {
                let mut atomic_contacts: Vec<ResultEntry> = vec![];

                // Clashes and VdW contacts
                let vdw =
                    find_vdw_contact(&x.0, &x.1, self.vdw_comp_factor).map(|intxn| ResultEntry {
                        interaction: intxn,
                        ligand: hierarchy_to_entity(&x.0),
                        receptor: hierarchy_to_entity(&x.1),
                        distance: x.0.atom().distance(x.1.atom()),
                    });
                atomic_contacts.extend(vdw.clone());

                // Skip checking for other interactions if there is a clash
                if vdw.is_some_and(|x| x.interaction == Interaction::StericClash) {
                    return Some(atomic_contacts);
                }

                // Hydrogen bonds and polar contacts
                let hbonds =
                    find_hydrogen_bond(&x.0, &x.1, self.vdw_comp_factor).map(|intxn| ResultEntry {
                        interaction: intxn,
                        ligand: hierarchy_to_entity(&x.0),
                        receptor: hierarchy_to_entity(&x.1),
                        distance: x.0.atom().distance(x.1.atom()),
                    });
                atomic_contacts.extend(hbonds);

                // C-H...O bonds
                let weak_hbonds =
                    find_weak_hydrogen_bond(&x.0, &x.1, self.vdw_comp_factor).map(|intxn| {
                        ResultEntry {
                            interaction: intxn,
                            ligand: hierarchy_to_entity(&x.0),
                            receptor: hierarchy_to_entity(&x.1),
                            distance: x.0.atom().distance(x.1.atom()),
                        }
                    });
                atomic_contacts.extend(weak_hbonds);

                // Ionic bonds
                let ionic_bonds = find_ionic_bond(&x.0, &x.1).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: hierarchy_to_entity(&x.0),
                    receptor: hierarchy_to_entity(&x.1),
                    distance: x.0.atom().distance(x.1.atom()),
                });
                atomic_contacts.extend(ionic_bonds);

                // Hydrophobic contacts
                let hydrophobic_contacts =
                    find_hydrophobic_contact(&x.0, &x.1).map(|intxn| ResultEntry {
                        interaction: intxn,
                        ligand: hierarchy_to_entity(&x.0),
                        receptor: hierarchy_to_entity(&x.1),
                        distance: x.0.atom().distance(x.1.atom()),
                    });
                atomic_contacts.extend(hydrophobic_contacts);

                Some(atomic_contacts)
            })
            .flatten()
            .collect::<Vec<ResultEntry>>()
    }

    fn get_ring_contacts(&self) {
        todo!()
    }
}

/// Find the absolute index of a residue in each chain.
///
/// Returns two mappings, one from residue to index, and one from index to residue.
fn build_residue_index(model: &PDB) -> HashMap<String, HashMap<(isize, String), usize>> {
    let mut res2idx: HashMap<String, HashMap<(isize, String), usize>> = HashMap::new();

    model.chains().for_each(|chain| {
        let chain_id = chain.id().to_string();
        res2idx.insert(
            chain_id.clone(),
            chain
                .residues()
                .enumerate()
                .map(|(i, residue)| {
                    let (resi, insertion) = residue.id();
                    let altloc = match insertion {
                        Some(insertion) => insertion.to_string(),
                        None => "".to_string(),
                    };
                    ((resi, altloc), i)
                })
                .collect::<HashMap<(isize, String), usize>>(),
        );
    });
    res2idx
}

fn build_ring_positions(
    model: &PDB,
) -> HashMap<String, HashMap<(isize, String), (na::Vector3<f64>, na::Vector3<f64>)>> {
    let mut rings: HashMap<String, HashMap<(isize, String), (na::Vector3<f64>, na::Vector3<f64>)>> =
        HashMap::new();

    model.chains().for_each(|chain| {
        let chain_id = chain.id().to_string();
        let ring_res = HashSet::from(["HIS", "PHE", "TYR", "TRP"]);
        rings.insert(
            chain_id.clone(),
            chain
                .residues()
                .filter(|r| ring_res.contains(r.name().unwrap()))
                .map(|residue| {
                    let (resi, insertion) = residue.id();
                    let altloc = match insertion {
                        Some(insertion) => insertion.to_string(),
                        None => "".to_string(),
                    };
                    ((resi, altloc), residue.ring_center_and_normal().unwrap())
                })
                .collect::<HashMap<(isize, String), (na::Vector3<f64>, na::Vector3<f64>)>>(),
        );
    });
    rings
}
