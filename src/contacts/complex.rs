use super::{
    InteractingEntity, Interaction, ProtonationMode, ResultEntry, find_cation_pi,
    find_hydrophobic_contact, find_ionic_bond_with_protonation,
    find_ionic_repulsion_with_protonation, find_pi_pi, find_vdw_contact,
    hbond::{find_hydrogen_bond, find_weak_hydrogen_bond},
    residues::{Plane, ResidueExt, ResidueId},
};
use crate::{
    StructureMetadata,
    metadata::{ExplicitBondKind, ResolvedBonds},
    structure::parse_groups,
};
use pdbtbx::*;
use rayon::prelude::*;
use rstar::{RTree, primitives::GeomWithData};
use std::collections::{HashMap, HashSet};

const MAX_PEPTIDE_C_N_DISTANCE: f64 = 1.8;
const GROUP1: u8 = 0b01;
const GROUP2: u8 = 0b10;

pub(crate) type SidechainStat<'a> = ((ResidueId<'a>, ResidueId<'a>), (f64, f64, f64));

struct IndexedAtom<'a> {
    entity: AtomConformerResidueChainModel<'a>,
    group: u8,
    residue_index: usize,
}

struct IndexedRing<'a> {
    residue: ResidueId<'a>,
    plane: Plane,
    group: u8,
    residue_index: usize,
}

/// The workhorse struct for identifying interactions in the model
pub struct InteractionComplex<'a> {
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

    atoms: Vec<IndexedAtom<'a>>,
    atom_tree: RTree<GeomWithData<[f64; 3], usize>>,
    /// Residue pairs connected by observed peptide-bond geometry.
    peptide_neighbors: Vec<(usize, usize)>,
    /// Maps ring residues to ring centers and normals
    rings: Vec<IndexedRing<'a>>,
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
        let peptide_neighbors = build_peptide_neighbors(model, metadata, &res2idx);
        let atoms = build_indexed_atoms(model, &ligand, &receptor, &res2idx);
        let atom_tree = RTree::bulk_load(
            atoms
                .iter()
                .enumerate()
                .map(|(index, atom)| {
                    let position = atom.entity.atom().pos();
                    GeomWithData::new([position.0, position.1, position.2], index)
                })
                .collect(),
        );
        let resolved_bonds = metadata.map_or_else(ResolvedBonds::default, |metadata| {
            metadata.resolve_bonds(model)
        });

        // Build a mapping of ring residue names to ring centers and normals
        let (rings, ring_err) = build_ring_positions(model, &ligand, &receptor, &res2idx);

        // Similarly, build a mapping of side chain planes
        let sc_planes = build_sc_plane_positions(model);

        Ok((
            Self {
                ligand,
                receptor,
                vdw_comp_factor,
                interacting_threshold,
                resolved_bonds,
                protonation,
                ph,
                atoms,
                atom_tree,
                peptide_neighbors,
                rings,
                sc_planes,
            },
            ring_err,
        ))
    }

    fn explicit_bond_kind(
        &self,
        first: &AtomConformerResidueChainModel,
        second: &AtomConformerResidueChainModel,
    ) -> Option<ExplicitBondKind> {
        self.resolved_bonds.kind_entities(first, second)
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
    fn should_compare_entities(&self, e1: &IndexedAtom, e2: &IndexedAtom, symmetric: bool) -> bool {
        // Ignore if any of the atoms is a hydrogen atom
        if (e1.entity.atom().element() == Some(&Element::H))
            | (e2.entity.atom().element() == Some(&Element::H))
        {
            return false;
        }

        let e1_res = ResidueId::from_hier(&e1.entity);
        let e2_res = ResidueId::from_hier(&e2.entity);
        self.should_compare_residues(
            &e1_res,
            &e2_res,
            e1.group,
            e2.group,
            e1.residue_index,
            e2.residue_index,
            symmetric,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn should_compare_residues(
        &self,
        r1: &ResidueId,
        r2: &ResidueId,
        group1: u8,
        group2: u8,
        residue1: usize,
        residue2: usize,
        symmetric: bool,
    ) -> bool {
        // Ignore if the two atoms are from different models
        if r1.model != r2.model {
            return false;
        }

        // Ignore if they are not a valid ligand-receptor pair
        if !(((group1 & GROUP1 != 0) && (group2 & GROUP2 != 0))
            | ((group2 & GROUP1 != 0) && (group1 & GROUP2 != 0)))
        {
            return false;
        }

        // Ignore if they are neighboring residues in the same chain
        if r1.chain == r2.chain {
            if residue1 == residue2
                || self
                    .peptide_neighbors
                    .binary_search(&ordered_index_pair(residue1, residue2))
                    .is_ok()
            {
                return false;
            }

            if symmetric { residue1 < residue2 } else { true }
        } else {
            // Across two chains, avoid duplicate comparisons when the chains exist on both sides,
            // e.g. H,A,B/H,A where H-A and A-H are the same interactions
            !(symmetric
                && group1 == GROUP1 | GROUP2
                && group2 == GROUP1 | GROUP2
                && (r1.chain > r2.chain))
        }
    }

    pub(crate) fn collect_sc_stats<'b>(
        &self,
        atomic_contacts: &'b [ResultEntry],
        ring_contacts: &'b [ResultEntry],
    ) -> Vec<SidechainStat<'b>> {
        let mut pairs = atomic_contacts
            .iter()
            .chain(ring_contacts)
            .filter_map(|contact| {
                let res1 = ResidueId::new(
                    contact.model,
                    contact.ligand.chain.as_str(),
                    contact.ligand.resi,
                    contact.ligand.insertion.as_str(),
                    contact.ligand.altloc.as_str(),
                    contact.ligand.resn.as_str(),
                );
                let res2 = ResidueId::new(
                    contact.model,
                    contact.receptor.chain.as_str(),
                    contact.receptor.resi,
                    contact.receptor.insertion.as_str(),
                    contact.receptor.altloc.as_str(),
                    contact.receptor.resn.as_str(),
                );
                (self.sc_planes.contains_key(&res1) && self.sc_planes.contains_key(&res2))
                    .then_some((res1, res2))
            })
            .collect::<Vec<_>>();
        pairs.sort_unstable();
        pairs.dedup();
        pairs
            .into_par_iter()
            .map(|(res1, res2)| {
                let res1_plane = &self.sc_planes[&res1];
                let res2_plane = &self.sc_planes[&res2];
                let centroid_dist = res1_plane.point_vec_dist(&res2_plane.center);
                let dihedral = res1_plane.dihedral(res2_plane);
                let centroid_angle = res1_plane.point_vec_angle(&res2_plane.center);
                ((res1, res2), (centroid_dist, dihedral, centroid_angle))
            })
            .collect()
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
        let max_radius_squared = self.interacting_threshold * self.interacting_threshold;

        // Find all atoms within the radius of the ligand atoms
        let ligand_neighbors: Vec<(usize, usize)> = self
            .atoms
            .iter()
            .enumerate()
            .filter(|(_, atom)| atom.group & GROUP1 != 0)
            .flat_map(|(first_index, first)| {
                self.atom_tree
                    .locate_within_distance(
                        {
                            let position = first.entity.atom().pos();
                            [position.0, position.1, position.2]
                        },
                        max_radius_squared,
                    )
                    .filter_map(|point| {
                        let second = &self.atoms[point.data];
                        (second.group & GROUP2 != 0
                            && self.should_compare_entities(first, second, true))
                        .then_some((first_index, point.data))
                    })
                    .collect::<Vec<_>>()
            })
            .collect();

        ligand_neighbors
            .par_iter()
            .map(|&(first_index, second_index)| {
                let e1 = &self.atoms[first_index].entity;
                let e2 = &self.atoms[second_index].entity;
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
                let vdw = find_vdw_contact(
                    e1,
                    e2,
                    self.vdw_comp_factor,
                    self.explicit_bond_kind(e1, e2),
                )
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
        let max_radius_squared = self.interacting_threshold * self.interacting_threshold;

        // Find ring - atom contacts
        let ring_atom_neighbors = self
            .rings
            .iter()
            .enumerate()
            .flat_map(|(ring_index, ring)| {
                self.atom_tree
                    .locate_within_distance(
                        [
                            ring.plane.center.x,
                            ring.plane.center.y,
                            ring.plane.center.z,
                        ],
                        max_radius_squared,
                    )
                    .filter_map(|point| {
                        let atom = &self.atoms[point.data];
                        let atom_residue = ResidueId::from_hier(&atom.entity);
                        self.should_compare_residues(
                            &ring.residue,
                            &atom_residue,
                            ring.group,
                            atom.group,
                            ring.residue_index,
                            atom.residue_index,
                            false,
                        )
                        .then_some((ring_index, point.data))
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        // Find ring-atom interactions
        ring_atom_neighbors
            .par_iter()
            .map(|&(ring_index, atom_index)| {
                let indexed_ring = &self.rings[ring_index];
                let k = &indexed_ring.residue;
                let ring = &indexed_ring.plane;
                let y = &self.atoms[atom_index].entity;
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
            .enumerate()
            .flat_map(|(first_index, first)| {
                self.rings
                    .iter()
                    .enumerate()
                    .filter_map(|(second_index, second)| {
                        self.should_compare_residues(
                            &first.residue,
                            &second.residue,
                            first.group,
                            second.group,
                            first.residue_index,
                            second.residue_index,
                            true,
                        )
                        .then_some((first_index, second_index))
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        // Find ring-ring interactions
        ring_ring_neighbors
            .par_iter()
            .filter_map(|&(first_index, second_index)| {
                let first = &self.rings[first_index];
                let second = &self.rings[second_index];
                let k1 = &first.residue;
                let k2 = &second.residue;
                let ring1 = &first.plane;
                let ring2 = &second.plane;
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
    let mut result = HashMap::new();
    let mut index = 0;
    for model in model.models() {
        for chain in model.chains() {
            for residue in chain.residues() {
                let (resi, insertion) = residue.id();
                for conformer in residue.conformers() {
                    result.insert(
                        ResidueId::new(
                            model.serial_number(),
                            chain.id(),
                            resi,
                            insertion.unwrap_or(""),
                            conformer.alternative_location().unwrap_or(""),
                            residue.name().unwrap_or(""),
                        ),
                        index,
                    );
                }
                index += 1;
            }
        }
    }
    result
}

fn group_mask(chain: &str, ligand: &HashSet<String>, receptor: &HashSet<String>) -> u8 {
    u8::from(ligand.contains(chain)) * GROUP1 + u8::from(receptor.contains(chain)) * GROUP2
}

fn build_indexed_atoms<'a>(
    model: &'a PDB,
    ligand: &HashSet<String>,
    receptor: &HashSet<String>,
    residues: &HashMap<ResidueId<'a>, usize>,
) -> Vec<IndexedAtom<'a>> {
    model
        .atoms_with_hierarchy()
        .filter(|entity| entity.atom().element() != Some(&Element::H))
        .filter_map(|entity| {
            let group = group_mask(entity.chain().id(), ligand, receptor);
            (group != 0).then(|| IndexedAtom {
                residue_index: residues[&ResidueId::from_hier(&entity)],
                entity,
                group,
            })
        })
        .collect()
}

fn build_peptide_neighbors<'a>(
    pdb: &'a PDB,
    metadata: Option<&StructureMetadata>,
    residue_indices: &HashMap<ResidueId<'a>, usize>,
) -> Vec<(usize, usize)> {
    let mut neighbors = Vec::new();
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
                    let residue_index = |residue: &Residue| {
                        let conformer = residue.conformers().next().unwrap();
                        residue_indices[&ResidueId::new(
                            model_id,
                            chain_id,
                            residue.serial_number(),
                            residue.insertion_code().unwrap_or(""),
                            conformer.alternative_location().unwrap_or(""),
                            residue.name().unwrap_or(""),
                        )]
                    };
                    let first_index = residue_index(first);
                    let second_index = residue_index(second);
                    neighbors.push(ordered_index_pair(first_index, second_index));
                }
            }
        }
    }
    neighbors.sort_unstable();
    neighbors.dedup();
    neighbors
}

fn ordered_index_pair(first: usize, second: usize) -> (usize, usize) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn build_ring_positions<'a>(
    model: &'a PDB,
    ligand: &HashSet<String>,
    receptor: &HashSet<String>,
    residues: &HashMap<ResidueId<'a>, usize>,
) -> (Vec<IndexedRing<'a>>, Vec<String>) {
    let mut ring_positions = Vec::new();
    let mut errors = Vec::new();

    for m in model.models() {
        let model_id = m.serial_number();
        for c in m.chains() {
            let chain_id = c.id();
            for r in c.residues().filter(|r| {
                r.name().is_some_and(|name| {
                    matches!(
                        name,
                        "HIS"
                            | "HID"
                            | "HIE"
                            | "HIP"
                            | "HSD"
                            | "HSE"
                            | "HSP"
                            | "PHE"
                            | "TYR"
                            | "TRP"
                    )
                })
            }) {
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
                            ring_positions.push(IndexedRing {
                                group: group_mask(chain_id, ligand, receptor),
                                residue_index: residues[&res_id],
                                residue: res_id,
                                plane: ring,
                            });
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

        let chains = HashSet::from(["A".to_string()]);
        let residues = build_residue_index(&pdb);
        let (rings, errors) = build_ring_positions(&pdb, &chains, &chains, &residues);

        assert!(errors.is_empty());
        assert_eq!(rings.len(), 2);
        let first = &rings
            .iter()
            .find(|ring| ring.residue.model == 1)
            .unwrap()
            .plane;
        let second = &rings
            .iter()
            .find(|ring| ring.residue.model == 2)
            .unwrap()
            .plane;
        assert!((second.center.x - first.center.x) > 5.0);
    }
}
