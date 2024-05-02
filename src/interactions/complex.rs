use super::hbond::*;
use super::structs::Interaction;
use super::structs::{InteractingEntity, ResultEntry};
use super::vdw::*;
use crate::utils::parse_groups;

use pdbtbx::*;
use rayon::prelude::*;
use std::collections::HashSet;
use tracing::debug;

pub struct InteractionComplex {
    pub model: PDB,
    pub ligand: HashSet<String>,
    pub receptor: HashSet<String>,
    pub vdw_comp_factor: f64,
    pub interacting_threshold: f64,
}

impl InteractionComplex {
    pub fn new(model: PDB, groups: &str, vdw_comp_factor: f64, interacting_threshold: f64) -> Self {
        let all_chains: HashSet<String> = model.par_chains().map(|c| c.id().to_string()).collect();
        let (ligand, receptor) = parse_groups(&all_chains, groups);

        debug!("Parsed ligand chains {ligand:?}; receptor chains {receptor:?}");
        Self {
            model,
            ligand,
            receptor,
            vdw_comp_factor,
            interacting_threshold,
        }
    }
}

pub trait Interactions {
    fn get_ppi(&self);
    fn get_atomic_contacts(&self) -> Vec<ResultEntry>;
    fn get_ring_contacts(&self);
}

impl Interactions for InteractionComplex {
    fn get_ppi(&self) {
        // self.get_atomic_contacts()
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

                Some(atomic_contacts)
            })
            .flatten()
            .collect::<Vec<ResultEntry>>()
    }

    fn get_ring_contacts(&self) {
        todo!()
    }
}

fn hierarchy_to_entity(hierarchy: &AtomConformerResidueChainModel<'_>) -> InteractingEntity {
    InteractingEntity {
        chain: hierarchy.chain().id().to_string(),
        resn: hierarchy.residue().name().unwrap().to_string(),
        resi: hierarchy.residue().serial_number(),
        altloc: hierarchy
            .conformer()
            .alternative_location()
            .unwrap_or("")
            .to_string(),
        atomn: hierarchy.atom().name().to_string(),
        atomi: hierarchy.atom().serial_number(),
    }
}
