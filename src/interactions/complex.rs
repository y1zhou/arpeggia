use crate::interactions::hbond::*;
use crate::interactions::structs::{InteractingEntity, Interaction, ResultEntry};
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
