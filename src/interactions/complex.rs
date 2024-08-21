use super::{
    find_cation_pi, find_hydrogen_bond, find_hydrophobic_contact, find_ionic_bond,
    find_ionic_repulsion, find_pi_pi, find_vdw_contact, find_weak_hydrogen_bond, point_ring_dist,
    InteractingEntity, Interaction, ResultEntry,
};
use crate::{
    residues::{ResidueExt, ResidueId, Ring},
    utils::parse_groups,
};
use pdbtbx::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

/// The workhorse struct for identifying interactions in the model
pub struct InteractionComplex<'a> {
    /// All information present in the atomic model
    pub model: &'a PDB,
    /// All ligand chains
    pub ligand: HashSet<String>,
    /// All receptor chains
    pub receptor: HashSet<String>,
    /// Compensation factor for VdW radii dependent interaction types
    pub vdw_comp_factor: f64,
    /// Distance cutoff when searching for neighboring atoms
    pub interacting_threshold: f64,

    /// Maps residue names to unique indices
    res2idx: HashMap<ResidueId<'a>, usize>,
    /// Maps ring residues to ring centers and normals
    rings: HashMap<ResidueId<'a>, Ring>,
}

impl<'a> InteractionComplex<'a> {
    pub fn new(
        model: &'a PDB,
        groups: &'a str,
        vdw_comp_factor: f64,
        interacting_threshold: f64,
    ) -> Self {
        // Parse all chains and input chain groups
        let all_chains: HashSet<String> = model.par_chains().map(|c| c.id().to_string()).collect();
        let (ligand, receptor) = parse_groups(&all_chains, groups);

        // Build a mapping of residue names to indices
        let res2idx = build_residue_index(model);

        // Build a mapping of ring residue names to ring centers and normals
        let rings = build_ring_positions(model);

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

    fn is_neighboring_hierarchy(
        &self,
        x: &AtomConformerResidueChainModel,
        y: &AtomConformerResidueChainModel,
    ) -> bool {
        let x_id = ResidueId::from_hier(x);
        let y_id = ResidueId::from_hier(y);
        self.is_neighboring_res_pair(&x_id, &y_id)
    }

    fn is_neighboring_res_pair(&self, x: &ResidueId, y: &ResidueId) -> bool {
        if x.chain != y.chain {
            return false;
        }
        let x_idx = self.res2idx[x];
        let y_idx = self.res2idx[y];

        match x_idx {
            0 => (y_idx == x_idx) | (y_idx == x_idx + 1),
            _ => (y_idx == x_idx - 1) | (y_idx == x_idx) | (y_idx == x_idx + 1),
        }
    }
}

pub trait Interactions {
    /// Get all atomic interactions between the ligand and receptor.
    fn get_atomic_contacts(&self) -> Vec<ResultEntry>;

    /// Get all ring-atom interactions between the ligand and receptor.
    fn get_ring_atom_contacts(&self) -> Vec<ResultEntry>;
    /// Get all ring-ring interactions between the ligand and receptor.
    fn get_ring_ring_contacts(&self) -> Vec<ResultEntry>;
}

impl<'a> Interactions for InteractionComplex<'a> {
    fn get_atomic_contacts(&self) -> Vec<ResultEntry> {
        let tree = self.model.create_hierarchy_rtree();
        let max_radius_squared = self.interacting_threshold * self.interacting_threshold;

        // Find all atoms within the radius of the ligand atoms
        let ligand_neighbors: Vec<(
            AtomConformerResidueChainModel,
            &AtomConformerResidueChainModel,
        )> = self
            .model
            .atoms_with_hierarchy()
            .filter(|x| {
                self.ligand.contains(x.chain().id()) & (x.atom().element().unwrap() != &Element::H)
            })
            .flat_map(|x| {
                tree.locate_within_distance(x.atom().pos(), max_radius_squared)
                    .filter(|y| {
                        self.receptor.contains(y.chain().id())
                            & (y.atom().element().unwrap() != &Element::H)
                            // Skip lower resi if both entities are on the same chain
                            & !((x.chain().id() == y.chain().id())
                                & (x.residue().serial_number() > y.residue().serial_number()))
                            // Skip neighboring residues
                            & !self.is_neighboring_hierarchy(&x, y)
                    })
                    .map(|y| (x.clone(), y))
                    .collect::<Vec<(
                        AtomConformerResidueChainModel,
                        &AtomConformerResidueChainModel,
                    )>>()
            })
            .collect();

        ligand_neighbors
            .par_iter()
            .filter_map(|(e1, e2)| {
                let mut atomic_contacts: Vec<ResultEntry> = vec![];

                // Clashes and VdW contacts
                let vdw = find_vdw_contact(e1, e2, self.vdw_comp_factor).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: InteractingEntity::from_hier(e1),
                    receptor: InteractingEntity::from_hier(e2),
                    distance: e1.atom().distance(e2.atom()),
                });
                atomic_contacts.extend(vdw.clone());

                // Skip checking for other interactions if there is a clash
                if vdw.is_some_and(|x| x.interaction == Interaction::StericClash) {
                    return Some(atomic_contacts);
                }

                // Hydrogen bonds and polar contacts
                let hbonds =
                    find_hydrogen_bond(e1, e2, self.vdw_comp_factor).map(|intxn| ResultEntry {
                        interaction: intxn,
                        ligand: InteractingEntity::from_hier(e1),
                        receptor: InteractingEntity::from_hier(e2),
                        distance: e1.atom().distance(e2.atom()),
                    });
                atomic_contacts.extend(hbonds);

                // C-H...O bonds
                let weak_hbonds =
                    find_weak_hydrogen_bond(e1, e2, self.vdw_comp_factor).map(|intxn| {
                        ResultEntry {
                            interaction: intxn,
                            ligand: InteractingEntity::from_hier(e1),
                            receptor: InteractingEntity::from_hier(e2),
                            distance: e1.atom().distance(e2.atom()),
                        }
                    });
                atomic_contacts.extend(weak_hbonds);

                // Ionic bonds
                let ionic_bonds = find_ionic_bond(e1, e2).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: InteractingEntity::from_hier(e1),
                    receptor: InteractingEntity::from_hier(e2),
                    distance: e1.atom().distance(e2.atom()),
                });
                atomic_contacts.extend(ionic_bonds);

                // Charge-charge repulsions
                let charge_repulsions = find_ionic_repulsion(e1, e2).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: InteractingEntity::from_hier(e1),
                    receptor: InteractingEntity::from_hier(e2),
                    distance: e1.atom().distance(e2.atom()),
                });
                atomic_contacts.extend(charge_repulsions);

                // Hydrophobic contacts
                let hydrophobic_contacts =
                    find_hydrophobic_contact(e1, e2).map(|intxn| ResultEntry {
                        interaction: intxn,
                        ligand: InteractingEntity::from_hier(e1),
                        receptor: InteractingEntity::from_hier(e2),
                        distance: e1.atom().distance(e2.atom()),
                    });
                atomic_contacts.extend(hydrophobic_contacts);

                Some(atomic_contacts)
            })
            .flatten()
            .collect::<Vec<ResultEntry>>()
    }

    fn get_ring_atom_contacts(&self) -> Vec<ResultEntry> {
        let tree = self.model.create_hierarchy_rtree();
        let max_radius_squared = self.interacting_threshold * self.interacting_threshold;

        // Find ring - atom contacts
        let ring_atom_neighbors = self
            .rings
            .iter()
            .flat_map(|(x, v)| {
                tree.locate_within_distance(
                    (v.center.x, v.center.y, v.center.z),
                    max_radius_squared,
                )
                // Ring on ligand or receptor, and atom on the other side
                .filter(|y| {
                    (self.ligand.contains(x.chain) & self.receptor.contains(y.chain().id()))
                        | (self.ligand.contains(y.chain().id()) & self.receptor.contains(x.chain))
                })
                .filter(|y_hier| {
                    let y = ResidueId::from_hier(y_hier);
                    (y_hier.atom().element().unwrap() != &Element::H)
                        // Skip lower resi if both entities are on the same chain
                        & !((x.chain == y.chain)
                            & (x.resi > y.resi))
                        // Skip neighboring residues
                        & !self.is_neighboring_res_pair(x, &y)
                })
                .map(|y| (x, v, y))
                .collect::<Vec<(&ResidueId, &Ring, &AtomConformerResidueChainModel)>>()
            })
            .collect::<Vec<(&ResidueId, &Ring, &AtomConformerResidueChainModel)>>();

        // Find ring-atom interactions
        ring_atom_neighbors
            .par_iter()
            .filter_map(|(k, ring, y)| {
                let mut ring_contacts = Vec::new();

                // Cation-pi interactions
                let dist = point_ring_dist(ring, &y.atom().pos());
                let cation_pi_contacts = find_cation_pi(ring, y).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: InteractingEntity::new(k.chain, k.resi, k.altloc, k.resn, "Ring", 0),
                    receptor: InteractingEntity::from_hier(y),
                    distance: dist,
                });
                ring_contacts.extend(cation_pi_contacts);

                Some(ring_contacts)
            })
            .flatten()
            .collect::<Vec<ResultEntry>>()
    }

    fn get_ring_ring_contacts(&self) -> Vec<ResultEntry> {
        // Find ring - ring contacts
        let ring_ring_neighbors = self
            .rings
            .iter()
            .filter(|(k, _)| self.ligand.contains(k.chain))
            .flat_map(|(k1, ring1)| {
                self.rings
                    .iter()
                    .filter(
                        |(k2, _)| {
                            self.receptor.contains(k2.chain)
                        & !((k1.chain == k2.chain) & (k1.resi > k2.resi)) // Skip lower resi if both entities are on the same chain
                        & !self.is_neighboring_res_pair(k1, k2)
                        }, // Skip neighboring residues
                    )
                    .map(|(k2, ring2)| (k1, ring1, k2, ring2))
                    .collect::<Vec<(&ResidueId, &Ring, &ResidueId, &Ring)>>()
            })
            .collect::<Vec<(&ResidueId, &Ring, &ResidueId, &Ring)>>();

        // Find ring-ring interactions
        ring_ring_neighbors
            .par_iter()
            .filter_map(|(k1, ring1, k2, ring2)| {
                let dist = (ring1.center - ring2.center).norm();
                let pi_pi_contacts = find_pi_pi(ring1, ring2).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: InteractingEntity::new(
                        k1.chain, k1.resi, k1.altloc, k1.resn, "Ring", 0,
                    ),
                    receptor: InteractingEntity::new(
                        k2.chain, k2.resi, k2.altloc, k2.resn, "Ring", 0,
                    ),
                    distance: dist,
                });

                Some(pi_pi_contacts)
            })
            .flatten()
            .collect::<Vec<ResultEntry>>()
    }
}

/// Find the absolute index of a residue in each chain.
///
/// Returns two mappings, one from residue to index, and one from index to residue.
fn build_residue_index(model: &PDB) -> HashMap<ResidueId, usize> {
    model
        .chains()
        .flat_map(|c| {
            c.residues().enumerate().map(move |(i, residue)| {
                let res_id = ResidueId::from_residue(residue, c.id());

                (res_id, i)
            })
        })
        .collect::<HashMap<ResidueId, usize>>()
}

fn build_ring_positions(model: &PDB) -> HashMap<ResidueId, Ring> {
    let ring_res = HashSet::from(["HIS", "PHE", "TYR", "TRP"]);
    model
        .chains()
        .flat_map(|c| {
            c.residues()
                .filter(|r| ring_res.contains(r.name().unwrap()))
                .map(|r| {
                    (
                        ResidueId::from_residue(r, c.id()),
                        match r.ring_center_and_normal() {
                            Some(ring) => ring,
                            None => panic!("Failed to calculate ring position for residue {:?}", r),
                        },
                    )
                })
        })
        .collect()
}
