use crate::interactions::hbond::*;
use crate::interactions::structs::{InteractingEntity, Interaction, ResultEntry};
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
                ligand = model
                    .chains()
                    .map(|c| c.id().to_string())
                    .filter(|c| !receptor.contains(c))
                    .collect();
            } else {
                receptor = model
                    .chains()
                    .map(|c| c.id().to_string())
                    .filter(|c| !ligand.contains(c))
                    .collect();
            }
        }
        if ligand.is_empty() && receptor.is_empty() {
            panic!("No ligand or receptor chains passed!")
        }

        debug!("Parsed ligands {ligand:?}; receptors {receptor:?}");
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

        // Find all atoms within the radius of the ligand
        let ligand_neighbors: Vec<ResultEntry> = self
            .model
            .atoms_with_hierarchy()
            .filter(|x| self.ligand.contains(x.chain().id()))
            .filter(|x| is_hydrogen_donor(x.conformer(), x.atom()))
            .flat_map(|x| {
                // TODO: x.conformer().alternative_location() may also be needed
                let neighbor_entities: Vec<InteractingEntity> = tree
                    .locate_within_distance(x.atom().pos(), max_radius_squared)
                    .filter(|y| {
                        (y.chain().id() != x.chain().id())
                            && (is_hydrogen_acceptor(y.conformer(), y.atom()))
                    })
                    .map(|y| hierarchy_to_entity(y))
                    .collect();

                match neighbor_entities.is_empty() {
                    false => {
                        let x_entity = hierarchy_to_entity(&x);

                        neighbor_entities
                            .iter()
                            .map(|y| ResultEntry {
                                interaction: Interaction::HydrogenBond,
                                ligand: x_entity.clone(),
                                receptor: y.clone(),
                            })
                            .collect()
                    }
                    true => vec![],
                }
            })
            .collect();

        ligand_neighbors
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
