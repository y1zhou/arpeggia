use super::{
    InteractingEntity, Interaction, ProtonationMode, ResultEntry, find_cation_pi,
    find_hydrophobic_contact, find_ionic_bond_with_protonation,
    find_ionic_repulsion_with_protonation, find_pi_pi, find_vdw_contact,
    hbond::{find_hydrogen_bond, find_weak_hydrogen_bond},
    residues::{Plane, ResidueExt, ResidueId},
};
use crate::{StructureMetadata, metadata::ResolvedBonds, structure::parse_groups};
use pdbtbx::*;
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};

const MAX_PEPTIDE_C_N_DISTANCE: f64 = 1.8;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct ResiduePosition<'a> {
    model: usize,
    chain: &'a str,
    residue: isize,
    insertion: &'a str,
}

impl<'a> ResiduePosition<'a> {
    fn from_id(residue: &'a ResidueId<'a>) -> Self {
        Self {
            model: residue.model,
            chain: residue.chain,
            residue: residue.resi,
            insertion: residue.insertion,
        }
    }

    fn from_residue(model: usize, chain: &'a str, residue: &'a Residue) -> Self {
        Self {
            model,
            chain,
            residue: residue.serial_number(),
            insertion: residue.insertion_code().unwrap_or(""),
        }
    }
}

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
    resolved_bonds: ResolvedBonds,
    protonation: ProtonationMode,
    ph: f64,

    /// Maps residue names to unique indices
    res2idx: HashMap<ResidueId<'a>, usize>,
    /// Residue pairs connected by observed peptide-bond geometry.
    peptide_neighbors: HashSet<(ResiduePosition<'a>, ResiduePosition<'a>)>,
    /// Maps ring residues to ring centers and normals
    rings: HashMap<ResidueId<'a>, Plane>,
    /// Map residues to side chain planes
    sc_planes: HashMap<ResidueId<'a>, Plane>,
}

impl<'a> InteractionComplex<'a> {
    /// Construct a contact complex with all scientific policies explicit.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_options(
        model: &'a PDB,
        metadata: Option<&'a StructureMetadata>,
        groups: &'a str,
        vdw_comp_factor: f64,
        interacting_threshold: f64,
        protonation: ProtonationMode,
        ph: f64,
    ) -> crate::ArpeggiaResult<(Self, Vec<String>)> {
        // Parse all chains and input chain groups
        let all_chains: HashSet<String> = model.par_chains().map(|c| c.id().to_string()).collect();
        let (ligand, receptor) = parse_groups(&all_chains, groups)?;

        // Build a mapping of residue names to indices
        let res2idx = build_residue_index(model);
        let peptide_neighbors = build_peptide_neighbors(model, metadata);
        let resolved_bonds = metadata.map_or_else(ResolvedBonds::default, |metadata| {
            metadata.resolve_bonds(model)
        });

        // Build a mapping of ring residue names to ring centers and normals
        let (rings, ring_err) = build_ring_positions(model);

        // Similarly, build a mapping of side chain planes
        let sc_planes = build_sc_plane_positions(model);

        Ok((
            Self {
                model,
                ligand,
                receptor,
                vdw_comp_factor,
                interacting_threshold,
                resolved_bonds,
                protonation,
                ph,
                res2idx,
                peptide_neighbors,
                rings,
                sc_planes,
            },
            ring_err,
        ))
    }

    fn has_explicit_bond(
        &self,
        first: &AtomConformerResidueChainModel,
        second: &AtomConformerResidueChainModel,
    ) -> bool {
        self.resolved_bonds.contains_entities(first, second)
    }

    pub(crate) fn resolved_bonds(&self) -> &ResolvedBonds {
        &self.resolved_bonds
    }

    /// Determine if two entities need to be checked for interactions or not.
    /// For interactions within the same chain, skip if e1 appears later in the chain than e2,
    /// or if e2 is in the same or neighboring residue as e1.
    /// Across chains, first see if the two chains both appear as ligands and receptors.
    /// In such cases, we avoid calculations where c1 > c2 if the interaction is symmetric.
    /// Currently, only ring-atom interactions are asymmetric.
    fn should_compare_entities(
        &self,
        e1: &AtomConformerResidueChainModel,
        e2: &AtomConformerResidueChainModel,
        symmetric: bool,
    ) -> bool {
        // Ignore if any of the atoms is a hydrogen atom
        if (e1.atom().element() == Some(&Element::H)) | (e2.atom().element() == Some(&Element::H)) {
            return false;
        }

        let e1_res = ResidueId::from_hier(e1);
        let e2_res = ResidueId::from_hier(e2);
        self.should_compare_residues(&e1_res, &e2_res, symmetric)
    }

    fn should_compare_residues(&self, r1: &ResidueId, r2: &ResidueId, symmetric: bool) -> bool {
        // Ignore if the two atoms are from different models
        if r1.model != r2.model {
            return false;
        }

        // Ignore if they are not a valid ligand-receptor pair
        if !((self.ligand.contains(r1.chain) && self.receptor.contains(r2.chain))
            | (self.ligand.contains(r2.chain) && self.receptor.contains(r1.chain)))
        {
            return false;
        }

        // Ignore if they are neighboring residues in the same chain
        if r1.chain == r2.chain {
            let r1_position = ResiduePosition::from_id(r1);
            let r2_position = ResiduePosition::from_id(r2);
            if r1_position == r2_position
                || self.peptide_neighbors.contains(&(r1_position, r2_position))
            {
                return false;
            }
            let e1_idx = self.res2idx[r1];
            let e2_idx = self.res2idx[r2];

            if symmetric { e1_idx < e2_idx } else { true }
        } else {
            // Across two chains, avoid duplicate comparisons when the chains exist on both sides,
            // e.g. H,A,B/H,A where H-A and A-H are the same interactions
            !(symmetric
                && self.receptor.contains(r1.chain)
                && self.receptor.contains(r2.chain)
                && self.ligand.contains(r1.chain)
                && self.ligand.contains(r2.chain)
                && (r1.chain > r2.chain))
        }
    }

    pub(crate) fn collect_sc_stats(
        &'_ self,
        contacts: &'a [ResultEntry],
    ) -> HashMap<(ResidueId<'_>, ResidueId<'_>), (f64, f64, f64)> {
        contacts
            .par_iter()
            .filter_map(|contact| {
                let res1 = ResidueId::new(
                    contact.model,
                    contact.ligand.chain.as_str(),
                    contact.ligand.resi,
                    contact.ligand.insertion.as_str(),
                    contact.ligand.altloc.as_str(),
                    contact.ligand.resn.as_str(),
                );
                let res1_plane = self.sc_planes.get(&res1)?;
                let res2 = ResidueId::new(
                    contact.model,
                    contact.receptor.chain.as_str(),
                    contact.receptor.resi,
                    contact.receptor.insertion.as_str(),
                    contact.receptor.altloc.as_str(),
                    contact.receptor.resn.as_str(),
                );
                let res2_plane = self.sc_planes.get(&res2)?;
                let centroid_dist = res1_plane.point_vec_dist(&res2_plane.center);
                let dihedral = res1_plane.dihedral(res2_plane);
                let centroid_angle = res1_plane.point_vec_angle(&res2_plane.center);
                Some(((res1, res2), (centroid_dist, dihedral, centroid_angle)))
            })
            .collect::<HashMap<(ResidueId, ResidueId), (f64, f64, f64)>>()
    }
}

/// Trait for calculating PPIs.
pub trait Interactions {
    /// Get all atomic interactions between the ligand and receptor.
    fn get_atomic_contacts(&self) -> Vec<ResultEntry>;

    /// Get all ring-atom interactions between the ligand and receptor.
    fn get_ring_atom_contacts(&self) -> Vec<ResultEntry>;
    /// Get all ring-ring interactions between the ligand and receptor.
    fn get_ring_ring_contacts(&self) -> Vec<ResultEntry>;
}

impl Interactions for InteractionComplex<'_> {
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
                self.ligand.contains(x.chain().id()) && (x.atom().element() != Some(&Element::H))
            })
            .flat_map(|x| {
                tree.locate_within_distance(x.atom().pos(), max_radius_squared)
                    .filter(|y| self.receptor.contains(y.chain().id()))
                    .filter(|y| self.should_compare_entities(&x, y, true))
                    .map(|y| (x.clone(), y))
                    .collect::<Vec<(
                        AtomConformerResidueChainModel,
                        &AtomConformerResidueChainModel,
                    )>>()
            })
            .collect();

        ligand_neighbors
            .par_iter()
            .map(|(e1, e2)| {
                let mut atomic_contacts: Vec<ResultEntry> = vec![];
                let model_id = e1.model().serial_number();
                let dist = e1.atom().distance(e2.atom());
                let make_entry = |interaction| ResultEntry {
                    model: model_id,
                    interaction,
                    ligand: InteractingEntity::from_hier(e1),
                    receptor: InteractingEntity::from_hier(e2),
                    distance: dist,
                };

                // Clashes and VdW contacts
                let vdw =
                    find_vdw_contact(e1, e2, self.vdw_comp_factor, self.has_explicit_bond(e1, e2))
                        .map(&make_entry);
                let steric_clash = vdw
                    .as_ref()
                    .is_some_and(|entry| entry.interaction == Interaction::StericClash);
                atomic_contacts.extend(vdw);

                // Skip checking for other interactions if there is a clash
                if steric_clash {
                    return atomic_contacts;
                }

                // Ionic bonds, Hydrogen bonds and polar contacts
                let ionic_bonds =
                    find_ionic_bond_with_protonation(e1, e2, self.protonation, self.ph);
                let hbonds = find_hydrogen_bond(e1, e2, self.vdw_comp_factor, &self.resolved_bonds);
                let electrostatic = match (ionic_bonds, hbonds) {
                    (Some(Interaction::IonicBond), Some(Interaction::HydrogenBond)) => {
                        Some(Interaction::SaltBridge)
                    }
                    // Potential protonation cannot establish a definite salt bridge. For
                    // other polar contacts, retain the stronger ionic classification.
                    (Some(ionic), Some(_)) => Some(ionic),
                    (Some(ionic), None) => Some(ionic),
                    (None, Some(hbond)) => Some(hbond),
                    _ => None,
                }
                .map(&make_entry);
                atomic_contacts.extend(electrostatic);

                // C-H...O bonds
                let weak_hbonds =
                    find_weak_hydrogen_bond(e1, e2, self.vdw_comp_factor, &self.resolved_bonds)
                        .map(&make_entry);
                atomic_contacts.extend(weak_hbonds);

                // Charge-charge repulsions
                let charge_repulsions =
                    find_ionic_repulsion_with_protonation(e1, e2, self.protonation, self.ph)
                        .map(&make_entry);
                atomic_contacts.extend(charge_repulsions);

                // Hydrophobic contacts
                let hydrophobic_contacts = find_hydrophobic_contact(e1, e2).map(&make_entry);
                atomic_contacts.extend(hydrophobic_contacts);

                atomic_contacts
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
                .filter(|y| {
                    let y_res = ResidueId::from_hier(y);
                    self.should_compare_residues(x, &y_res, false)
                })
                .map(|y| (x, v, y))
                .collect::<Vec<(&ResidueId, &Plane, &AtomConformerResidueChainModel)>>()
            })
            .collect::<Vec<(&ResidueId, &Plane, &AtomConformerResidueChainModel)>>();

        // Find ring-atom interactions
        ring_atom_neighbors
            .par_iter()
            .map(|(k, ring, y)| {
                let mut ring_contacts = Vec::new();

                // Cation-pi interactions
                let dist = ring.point_dist(&y.atom().pos());
                let cation_pi_contacts =
                    find_cation_pi(ring, y, self.protonation, self.ph).map(|intxn| ResultEntry {
                        model: k.model,
                        interaction: intxn,
                        ligand: InteractingEntity::new(
                            k.chain,
                            k.resi,
                            k.insertion,
                            k.altloc,
                            k.resn,
                            "Ring",
                            0,
                        ),
                        receptor: InteractingEntity::from_hier(y),
                        distance: dist,
                    });
                ring_contacts.extend(cation_pi_contacts);

                ring_contacts
            })
            .flatten()
            .collect::<Vec<ResultEntry>>()
    }

    fn get_ring_ring_contacts(&self) -> Vec<ResultEntry> {
        // Find ring - ring contacts
        let ring_ring_neighbors = self
            .rings
            .iter()
            .flat_map(|(k1, ring1)| {
                self.rings
                    .iter()
                    .filter(|(k2, _)| {
                        self.ligand.contains(k1.chain)
                            && self.receptor.contains(k2.chain)
                            && self.should_compare_residues(k1, k2, true)
                    })
                    .map(|(k2, ring2)| (k1, ring1, k2, ring2))
                    .collect::<Vec<(&ResidueId, &Plane, &ResidueId, &Plane)>>()
            })
            .collect::<Vec<(&ResidueId, &Plane, &ResidueId, &Plane)>>();

        // Find ring-ring interactions
        ring_ring_neighbors
            .par_iter()
            .filter_map(|(k1, ring1, k2, ring2)| {
                let dist = (ring1.center - ring2.center).norm();
                let pi_pi_contacts = find_pi_pi(ring1, ring2).map(|intxn| ResultEntry {
                    model: k1.model,
                    interaction: intxn,
                    ligand: InteractingEntity::new(
                        k1.chain,
                        k1.resi,
                        k1.insertion,
                        k1.altloc,
                        k1.resn,
                        "Ring",
                        0,
                    ),
                    receptor: InteractingEntity::new(
                        k2.chain,
                        k2.resi,
                        k2.insertion,
                        k2.altloc,
                        k2.resn,
                        "Ring",
                        0,
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
/// Returns a mapping from residue to index.
fn build_residue_index(model: &'_ PDB) -> HashMap<ResidueId<'_>, usize> {
    model
        .models()
        .flat_map(|m| {
            let model_id = m.serial_number();
            m.chains().flat_map(move |c| {
                let chain_id = c.id();
                c.residues().enumerate().flat_map(move |(i, residue)| {
                    // All conformers in the same residue should share the same index
                    let (resi, insertion) = residue.id();
                    let insertion = insertion.unwrap_or("");
                    let resn = residue.name().unwrap_or("");

                    residue.conformers().map(move |conformer| {
                        let res_id = ResidueId::new(
                            model_id,
                            chain_id,
                            resi,
                            insertion,
                            conformer.alternative_location().unwrap_or(""),
                            resn,
                        );

                        (res_id, i)
                    })
                })
            })
        })
        .collect::<HashMap<ResidueId, usize>>()
}

fn build_peptide_neighbors<'a>(
    pdb: &'a PDB,
    metadata: Option<&StructureMetadata>,
) -> HashSet<(ResiduePosition<'a>, ResiduePosition<'a>)> {
    let mut neighbors = HashSet::new();
    for model in pdb.models() {
        let model_id = model.serial_number();
        for chain in model.chains() {
            let chain_id = chain.id();
            let residues = chain.residues().collect::<Vec<_>>();
            for pair in residues.windows(2) {
                let first = pair[0];
                let second = pair[1];
                let explicitly_broken = metadata.is_some_and(|metadata| {
                    metadata.has_chain_break_after(
                        model_id,
                        chain_id,
                        first.serial_number(),
                        first.insertion_code().unwrap_or(""),
                        first.name().unwrap_or(""),
                    )
                });
                let peptide_geometry = first
                    .atoms()
                    .find(|atom| atom.name() == "C")
                    .zip(second.atoms().find(|atom| atom.name() == "N"))
                    .is_some_and(|(carbon, nitrogen)| {
                        carbon.distance(nitrogen) <= MAX_PEPTIDE_C_N_DISTANCE
                    });
                if peptide_geometry && !explicitly_broken {
                    let first = ResiduePosition::from_residue(model_id, chain_id, first);
                    let second = ResiduePosition::from_residue(model_id, chain_id, second);
                    neighbors.insert((first, second));
                    neighbors.insert((second, first));
                }
            }
        }
    }
    neighbors
}

fn build_ring_positions(model: &'_ PDB) -> (HashMap<ResidueId<'_>, Plane>, Vec<String>) {
    let ring_res = HashSet::from([
        "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "PHE", "TYR", "TRP",
    ]);
    let mut ring_positions = HashMap::new();
    let mut errors = Vec::new();

    for m in model.models() {
        let model_id = m.serial_number();
        for c in m.chains() {
            let chain_id = c.id();
            for r in c
                .residues()
                .filter(|r| r.name().is_some_and(|name| ring_res.contains(name)))
            {
                let (resi, insertion_code) = r.id();
                let insertion_code = insertion_code.unwrap_or("");
                let resn = r.name().unwrap_or("");
                for conformer in r.conformers() {
                    let res_id = ResidueId::new(
                        model_id,
                        chain_id,
                        resi,
                        insertion_code,
                        conformer.alternative_location().unwrap_or(""),
                        resn,
                    );
                    match r.center_and_normal(Some(r.ring_atoms())) {
                        Some(ring) => {
                            ring_positions.insert(res_id, ring);
                        }
                        None => {
                            errors
                                .push(format!("Failed to calculate ring position for {res_id:?}"));
                        }
                    }
                }
            }
        }
    }
    (ring_positions, errors)
}

fn build_sc_plane_positions(model: &PDB) -> HashMap<ResidueId<'_>, Plane> {
    let mut sc_plane_positions = HashMap::new();

    for m in model.models() {
        let model_id = m.serial_number();
        for c in m.chains() {
            let chain_id = c.id();
            for r in c.residues() {
                let (resi, insertion_code) = r.id();
                let insertion_code = insertion_code.unwrap_or("");
                let resn = r.name().unwrap_or("");
                for conformer in r.conformers() {
                    let res_id = ResidueId::new(
                        model_id,
                        chain_id,
                        resi,
                        insertion_code,
                        conformer.alternative_location().unwrap_or(""),
                        resn,
                    );
                    if let Some(plane) = r.center_and_normal(Some(r.sc_plane_atoms())) {
                        sc_plane_positions.insert(res_id, plane);
                    }
                }
            }
        }
    }
    sc_plane_positions
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;

    #[test]
    fn ring_planes_are_model_local() {
        let input = b"MODEL        1\n\
ATOM      1  CG  PHE A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CD1 PHE A   1       1.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      3  CE1 PHE A   1       1.500   1.000   0.000  1.00 20.00           C  \n\
ATOM      4  CZ  PHE A   1       1.000   2.000   0.000  1.00 20.00           C  \n\
ATOM      5  CE2 PHE A   1       0.000   2.000   0.000  1.00 20.00           C  \n\
ATOM      6  CD2 PHE A   1      -0.500   1.000   0.000  1.00 20.00           C  \n\
ENDMDL\n\
MODEL        2\n\
ATOM      1  CG  PHE A   1      10.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      2  CD1 PHE A   1      11.000   0.000   0.000  1.00 20.00           C  \n\
ATOM      3  CE1 PHE A   1      11.500   1.000   0.000  1.00 20.00           C  \n\
ATOM      4  CZ  PHE A   1      11.000   2.000   0.000  1.00 20.00           C  \n\
ATOM      5  CE2 PHE A   1      10.000   2.000   0.000  1.00 20.00           C  \n\
ATOM      6  CD2 PHE A   1       9.500   1.000   0.000  1.00 20.00           C  \n\
ENDMDL\nEND\n";
        let (pdb, _) = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap();

        let (rings, errors) = build_ring_positions(&pdb);

        assert!(errors.is_empty());
        assert_eq!(rings.len(), 2);
        let first = rings.iter().find(|(id, _)| id.model == 1).unwrap().1;
        let second = rings.iter().find(|(id, _)| id.model == 2).unwrap().1;
        assert!((second.center.x - first.center.x) > 5.0);
    }
}
