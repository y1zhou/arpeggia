use pdbtbx::*;
// use rstar::RTree;
use std::collections::HashSet;

use crate::residues::{is_hydrogen_acceptor, is_hydrogen_donor};

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

        println!("Parsed ligands {ligand:?}; receptors {receptor:?}");
        Self {
            model,
            ligand,
            receptor,
            vdw_comp_factor,
            interacting_threshold,
        }
    }
}

#[allow(dead_code)]
enum Interaction {
    StericClash,
    Covalent,
    VanDerWaals,

    Electrostatic,
    Aromatic,

    Hydrophobic,
}

#[allow(dead_code)]
enum Electrostatic {
    Ionic,
    HydrogenBond,
    WeakHydrogenBond, // C-H...O hydrogen bond
    // HalogenBond
    Polar, // hydrogen bonding without angle terms
}

#[allow(dead_code)]
enum Aromatic {
    PiDisplaced, // staggered stacking, parallel displaced
    PiT,         // perpendicular T-shaped
    PiSandwich,  // direct stacking, repulsive
    CationPi,
}

pub trait Interactions {
    fn get_ppi(&self);
    fn get_atomic_contacts(&self) -> HashSet<usize>;
    fn get_ring_contacts(&self);
}

impl Interactions for InteractionComplex {
    fn get_ppi(&self) {
        // self.get_atomic_contacts()
    }

    fn get_atomic_contacts(&self) -> HashSet<usize> {
        let tree = self.model.create_hierarchy_rtree();

        // Find all atoms within the radius of the ligand
        let ligand_neighbors: HashSet<usize> = self
            .model
            .atoms_with_hierarchy()
            .filter(|x| self.ligand.contains(x.chain().id()))
            .filter(|x| is_hydrogen_donor(x.residue(), x.atom()))
            .flat_map(|x| {
                // TODO: x.conformer().alternative_location() may also be needed
                let neighbor_atoms_idx: HashSet<usize> = tree
                    .locate_within_distance(
                        x.atom().pos(),
                        self.interacting_threshold * self.interacting_threshold,
                    )
                    .filter(|y| {
                        (y.chain().id() != x.chain().id())
                            && (is_hydrogen_acceptor(y.residue(), y.atom()))
                    })
                    .map(|y| {
                        println!(
                            "Chain {}, {}{} Atom {} ({})",
                            y.chain().id(),
                            y.residue().name().unwrap(),
                            y.residue().serial_number(),
                            y.atom().serial_number(),
                            y.atom().name()
                        );
                        y.atom().serial_number()
                    })
                    .collect();

                if !neighbor_atoms_idx.is_empty() {
                    println!(
                        "└─Chain {}, {}{} Atom {} ({}) neighbors",
                        x.chain().id(),
                        x.residue().name().unwrap(),
                        x.residue().serial_number(),
                        x.atom().serial_number(),
                        x.atom().name()
                    );
                }

                neighbor_atoms_idx
            })
            .collect();

        ligand_neighbors
    }

    fn get_ring_contacts(&self) {
        todo!()
    }
}
