use super::{
    find_cation_pi, find_hydrogen_bond, find_hydrophobic_contact, find_ionic_bond, find_pi_pi,
    find_vdw_contact, find_weak_hydrogen_bond, point_ring_dist, InteractingEntity, Interaction,
    ResultEntry,
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
    /// Maps ring residues to ring centers and normals
    rings: HashMap<(String, isize, String, String), Ring>,
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

    fn is_neighboring_hierarchy(
        &self,
        x: &AtomConformerResidueChainModel,
        y: &AtomConformerResidueChainModel,
    ) -> bool {
        let (x_resi, x_insertion) = x.residue().id();
        let x_altloc = x_insertion.unwrap_or("");
        let (y_resi, y_insertion) = y.residue().id();
        let y_altloc = y_insertion.unwrap_or("");
        self.is_neighboring_res_pair(
            x.chain().id(),
            x_resi,
            x_altloc,
            y.chain().id(),
            y_resi,
            y_altloc,
        )
    }

    fn is_neighboring_res_pair(
        &self,
        x_chain: &str,
        x_resi: isize,
        x_altloc: &str,
        y_chain: &str,
        y_resi: isize,
        y_altloc: &str,
    ) -> bool {
        if x_chain != y_chain {
            return false;
        }
        let x_idx = self.res2idx[&x_chain.to_string()][&(x_resi, x_altloc.to_string())];
        let y_idx = self.res2idx[&y_chain.to_string()][&(y_resi, y_altloc.to_string())];

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

impl Interactions for InteractionComplex {
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
                    ligand: hierarchy_to_entity(e1),
                    receptor: hierarchy_to_entity(e2),
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
                        ligand: hierarchy_to_entity(e1),
                        receptor: hierarchy_to_entity(e2),
                        distance: e1.atom().distance(e2.atom()),
                    });
                atomic_contacts.extend(hbonds);

                // C-H...O bonds
                let weak_hbonds =
                    find_weak_hydrogen_bond(e1, e2, self.vdw_comp_factor).map(|intxn| {
                        ResultEntry {
                            interaction: intxn,
                            ligand: hierarchy_to_entity(e1),
                            receptor: hierarchy_to_entity(e2),
                            distance: e1.atom().distance(e2.atom()),
                        }
                    });
                atomic_contacts.extend(weak_hbonds);

                // Ionic bonds
                let ionic_bonds = find_ionic_bond(e1, e2).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: hierarchy_to_entity(e1),
                    receptor: hierarchy_to_entity(e2),
                    distance: e1.atom().distance(e2.atom()),
                });
                atomic_contacts.extend(ionic_bonds);

                // Hydrophobic contacts
                let hydrophobic_contacts =
                    find_hydrophobic_contact(e1, e2).map(|intxn| ResultEntry {
                        interaction: intxn,
                        ligand: hierarchy_to_entity(e1),
                        receptor: hierarchy_to_entity(e2),
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
            .flat_map(|(k, v)| {
                tree
                    .locate_within_distance((v.center.x, v.center.y, v.center.z), max_radius_squared)
                    // Ring on ligand or receptor, and atom on the other side
                    .filter(|y|
                        (self.ligand.contains(k.0.as_str())
                            & self.receptor.contains(y.chain().id())
                        ) |
                        (self.ligand.contains(y.chain().id())
                            & self.receptor.contains(k.0.as_str())
                        )
                    )
                    .filter(|y|
                        (y.atom().element().unwrap() != &Element::H)
                        // Skip lower resi if both entities are on the same chain
                        & !((k.0 == y.chain().id())
                            & (k.1 > y.residue().serial_number()))
                        // Skip neighboring residues
                        & !self.is_neighboring_res_pair(k.0.as_str(), k.1, k.2.as_str(), y.chain().id(), y.residue().serial_number(), y.residue().insertion_code().unwrap_or(""))
                    ).map(|y| (k, v, y))
                    .collect::<Vec<(&(String, isize, String, String), &Ring, &AtomConformerResidueChainModel)>>()
            }).collect::<Vec<(&(String, isize, String, String), &Ring, &AtomConformerResidueChainModel)>>();

        // Find ring-atom interactions
        ring_atom_neighbors
            .par_iter()
            .filter_map(|(k, ring, x)| {
                let mut ring_contacts = Vec::new();

                // Cation-pi interactions
                let dist = point_ring_dist(ring, &x.atom().pos());
                let cation_pi_contacts = find_cation_pi(ring, x).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: InteractingEntity {
                        chain: k.0.clone(),
                        resi: k.1,
                        altloc: k.2.clone(),
                        resn: k.3.clone(),
                        atomn: "Ring".to_string(),
                        atomi: 0,
                    },
                    receptor: hierarchy_to_entity(x),
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
            .filter(|(k, _)| self.ligand.contains(k.0.as_str()))
            .flat_map(|(k1, ring1)| {
                self.rings.iter()
                    .filter(|(k2, _)|
                        self.receptor.contains(k2.0.as_str())
                        & !((k1.0 == k2.0) & (k1.1 > k2.1)) // Skip lower resi if both entities are on the same chain
                        & !self.is_neighboring_res_pair(k1.0.as_str(), k1.1, k1.2.as_str(), k2.0.as_str(), k2.1, k2.2.as_str()) // Skip neighboring residues
                    ).map(|(k2, ring2)| (k1, ring1, k2, ring2))
                    .collect::<Vec<(&(String, isize, String, String), &Ring, &(String, isize, String, String), &Ring)>>()
            }).collect::<Vec<(&(String, isize, String, String), &Ring, &(String, isize, String, String), &Ring)>>();

        // Find ring-ring interactions
        ring_ring_neighbors
            .par_iter()
            .filter_map(|(k1, ring1, k2, ring2)| {
                let dist = (ring1.center - ring2.center).norm();
                let pi_pi_contacts = find_pi_pi(ring1, ring2).map(|intxn| ResultEntry {
                    interaction: intxn,
                    ligand: InteractingEntity {
                        chain: k1.0.clone(),
                        resi: k1.1,
                        altloc: k1.2.clone(),
                        resn: k1.3.clone(),
                        atomn: "Ring".to_string(),
                        atomi: 0,
                    },
                    receptor: InteractingEntity {
                        chain: k2.0.clone(),
                        resi: k2.1,
                        altloc: k2.2.clone(),
                        resn: k2.3.clone(),
                        atomn: "Ring".to_string(),
                        atomi: 0,
                    },
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

fn build_ring_positions(model: &PDB) -> HashMap<(String, isize, String, String), Ring> {
    let ring_res = HashSet::from(["HIS", "PHE", "TYR", "TRP"]);
    model
        .atoms_with_hierarchy()
        .filter(|x| ring_res.contains(x.residue().name().unwrap()))
        .map(|x| {
            let (resi, insertion) = x.residue().id();
            let altloc = match insertion {
                Some(insertion) => insertion.to_string(),
                None => "".to_string(),
            };

            (
                (
                    x.chain().id().to_string(),
                    resi,
                    altloc,
                    x.residue().name().unwrap().to_string(),
                ),
                x.residue().ring_center_and_normal().unwrap(),
            )
        })
        .collect::<HashMap<(String, isize, String, String), Ring>>()
}
