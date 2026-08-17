use super::structs::Interaction;
use crate::{BondEndpoint, StructureMetadata};

use pdbtbx::*;
use rayon::prelude::*;
use std::collections::HashSet;

const HYDROGEN_BOND_DIST: f64 = 4.0;
const POLAR_DIST: f64 = 3.5;

/// Search for hydrogen bonds and polar contacts.
///
/// ## Details
///
/// It first checks if the two input residues are hydrogen bond donors and acceptors.
/// If so, it collects hydrogens attached to that donor and checks if any satisfy:
///
/// * dist(H, acceptor) <= vdw_radius(H) + vdw_radius(acceptor) + vdw_comp_factor, and
/// * angle(donor, H, acceptor) >= 90°,
///
/// where `vdw_comp_factor` is the compensation factor for VdW radii dependent interactions.
/// If the conditions are met, it returns [`Interaction::HydrogenBond`].
/// If not, it checks for [`Interaction::PolarContact`] by: dist(donor, acceptor) <= [`POLAR_DIST`].
///
/// For further details, see:
///
/// * [PyMOL wiki](https://pymolwiki.org/index.php/Displaying_Biochemical_Properties#Hydrogen_bonds_and_Polar_Contacts)
/// * [Arpeggio implementation](https://github.com/PDBeurope/arpeggio/blob/258855b8ba13447f2776b232ca32884d637c6a9c/arpeggio/core/utils.py#L73)
///
/// [WARNING] Water-mediated hydrogen-bond geometry is not currently modeled.
/// Search strong hydrogen bonds using canonical names or explicit input bonds.
pub fn find_hydrogen_bond(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
    metadata: Option<&StructureMetadata>,
) -> Option<Interaction> {
    if let Some((donor, acceptor)) = is_donor_acceptor_pair(entity1, entity2) {
        let da_dist = donor.atom().distance(acceptor.atom());
        if da_dist <= HYDROGEN_BOND_DIST {
            let donor_h: Vec<&Atom> = donor
                .residue()
                .par_atoms()
                .filter(|atom| {
                    atom.element() == Some(&Element::H)
                        && (is_hydrogen_for_donor(
                            donor.conformer().name(),
                            donor.atom().name(),
                            atom.name(),
                        ) || metadata.is_some_and(|metadata| {
                            metadata.has_bond(
                                &endpoint_for_entity(donor),
                                &BondEndpoint::from_parts(
                                    donor.chain().id(),
                                    donor.residue().serial_number(),
                                    donor.residue().insertion_code().unwrap_or(""),
                                    atom.name(),
                                ),
                            )
                        }))
                })
                .collect();

            // Hydrogen bonds are stricter as they have angle restrictions
            let acceptor_vdw = acceptor.atom().element()?.atomic_radius().van_der_waals?;
            let h_vdw: f64 = Element::H.atomic_radius().van_der_waals.unwrap();
            if donor_h.par_iter().any(|h| {
                (h.distance(acceptor.atom()) <= h_vdw + acceptor_vdw + vdw_comp_factor)
                    && (donor.atom().angle(h, acceptor.atom()) >= 90.0)
            }) {
                return Some(Interaction::HydrogenBond);
            }
        }
        if da_dist <= POLAR_DIST {
            // Polar interactions are more relaxed as only distance is checked
            return Some(Interaction::PolarContact);
        }
    }
    None
}

/// Search for weak hydrogen bonds and weak polar contacts.
///
/// Overall the same as [`find_hydrogen_bond`], except for:
///
/// * Seeks for C-H...O bonds instead of the usual N and O donors.
/// * Checks for angle(donor, H, acceptor) >= 130° instead of 90°.
///
/// Search weak hydrogen bonds using canonical names or explicit input bonds.
pub fn find_weak_hydrogen_bond(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
    metadata: Option<&StructureMetadata>,
) -> Option<Interaction> {
    if let Some((donor, acceptor)) = is_weak_donor_acceptor_pair(entity1, entity2) {
        let da_dist = donor.atom().distance(acceptor.atom());
        if da_dist <= HYDROGEN_BOND_DIST {
            let donor_h: Vec<&Atom> = donor
                .residue()
                .par_atoms()
                .filter(|atom| {
                    atom.element() == Some(&Element::H)
                        && (is_hydrogen_for_carbon(donor.atom().name(), atom.name())
                            || metadata.is_some_and(|metadata| {
                                metadata.has_bond(
                                    &endpoint_for_entity(donor),
                                    &BondEndpoint::from_parts(
                                        donor.chain().id(),
                                        donor.residue().serial_number(),
                                        donor.residue().insertion_code().unwrap_or(""),
                                        atom.name(),
                                    ),
                                )
                            }))
                })
                .collect();

            // Hydrogen bonds are stricter as they have angle restrictions
            let acceptor_vdw = acceptor.atom().element()?.atomic_radius().van_der_waals?;
            let h_vdw: f64 = Element::H.atomic_radius().van_der_waals.unwrap();
            if donor_h.par_iter().any(|h| {
                (h.distance(acceptor.atom()) <= h_vdw + acceptor_vdw + vdw_comp_factor)
                    && (donor.atom().angle(h, acceptor.atom()) >= 130.0)
            }) {
                return Some(Interaction::WeakHydrogenBond);
            }
        }
        if da_dist <= POLAR_DIST {
            // Polar interactions are more relaxed as only distance is checked
            return Some(Interaction::WeakPolarContact);
        }
    }
    None
}

fn endpoint_for_entity(entity: &AtomConformerResidueChainModel) -> BondEndpoint {
    BondEndpoint::from_parts(
        entity.chain().id(),
        entity.residue().serial_number(),
        entity.residue().insertion_code().unwrap_or(""),
        entity.atom().name(),
    )
}

/// Determine if the two entities are a valid hydrogen bond donor-acceptor pair
fn is_donor_acceptor_pair<'a>(
    entity1: &'a AtomConformerResidueChainModel<'a>,
    entity2: &'a AtomConformerResidueChainModel<'a>,
) -> Option<(
    &'a AtomConformerResidueChainModel<'a>,
    &'a AtomConformerResidueChainModel<'a>,
)> {
    let e1_conformer = entity1.conformer().name();
    let e2_conformer = entity2.conformer().name();
    let e1_atom = entity1.atom().name();
    let e2_atom = entity2.atom().name();

    if is_hydrogen_donor(e1_conformer, e1_atom) && is_hydrogen_acceptor(e2_conformer, e2_atom) {
        Some((entity1, entity2))
    } else if is_hydrogen_donor(e2_conformer, e2_atom)
        && is_hydrogen_acceptor(e1_conformer, e1_atom)
    {
        Some((entity2, entity1))
    } else {
        None
    }
}

/// Determine if the atom in the residue is a hydrogen acceptor
fn is_hydrogen_acceptor(res_name: &str, atom_name: &str) -> bool {
    // all the carbonyl oxygens in the main chain and on the terminals
    let oxygens = HashSet::from(["O", "OXT"]);
    if oxygens.contains(atom_name) && res_name != "HOH" {
        return true;
    }
    matches!(
        (res_name, atom_name),
        ("ASN", "OD1")
            // | ("ASN", "ND2")
            | ("ASP", "OD1" | "OD2")
            | ("GLN", "OE1")
            | ("GLU", "OE1" | "OE2")
            | (
                "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP",
                "ND1" | "NE2"
            )
            | ("SER", "OG")
            | ("THR", "OG1")
            | ("TYR", "OH")
            | ("MET", "SD") // 10.1021/jz300207k and 10.1002/prot.22327
            | ("CYS", "SG") // 10.1002/prot.22327
    )
}

/// Determine if the atom in the residue is a hydrogen donor
fn is_hydrogen_donor(res_name: &str, atom_name: &str) -> bool {
    // All amide niteogens in the main chain except proline
    if atom_name == "N" && res_name != "PRO" {
        return true;
    }
    matches!(
        (res_name, atom_name),
        ("ARG", "NE" | "NH1" | "NH2")
            | ("ASN", "ND2")
            | ("GLN", "NE2")
            | (
                "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP",
                "ND1" | "NE2"
            )
            | ("LYS", "NZ")
            | ("SER", "OG")
            | ("THR", "OG1")
            | ("TRP", "NE1")
            | ("TYR", "OH")
            | ("CYS", "SG") // 10.1002/prot.22327
    )
}

fn is_hydrogen_for_donor(res_name: &str, donor_name: &str, hydrogen_name: &str) -> bool {
    if donor_name == "N" {
        return matches!(
            hydrogen_name,
            "H" | "H1" | "H2" | "H3" | "HN" | "1H" | "2H" | "3H"
        );
    }
    match (res_name, donor_name) {
        ("ARG", "NE") => hydrogen_name == "HE",
        ("ARG", "NH1") => matches!(hydrogen_name, "HH11" | "HH12"),
        ("ARG", "NH2") => matches!(hydrogen_name, "HH21" | "HH22"),
        ("ASN", "ND2") => matches!(hydrogen_name, "HD21" | "HD22"),
        ("GLN", "NE2") => matches!(hydrogen_name, "HE21" | "HE22"),
        ("HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP", "ND1") => hydrogen_name == "HD1",
        ("HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP", "NE2") => hydrogen_name == "HE2",
        ("LYS", "NZ") => matches!(hydrogen_name, "HZ1" | "HZ2" | "HZ3"),
        ("SER", "OG") | ("CYS", "SG") => hydrogen_name == "HG",
        ("THR", "OG1") => hydrogen_name == "HG1",
        ("TRP", "NE1") => hydrogen_name == "HE1",
        ("TYR", "OH") => hydrogen_name == "HH",
        _ => false,
    }
}

pub(crate) fn count_donors_without_explicit_hydrogen(
    pdb: &PDB,
    metadata: Option<&StructureMetadata>,
) -> usize {
    pdb.atoms_with_hierarchy()
        .filter(|entity| {
            let residue_name = entity.conformer().name();
            let donor_name = entity.atom().name();
            is_hydrogen_donor(residue_name, donor_name)
                && !entity.residue().atoms().any(|atom| {
                    atom.element() == Some(&Element::H)
                        && (is_hydrogen_for_donor(residue_name, donor_name, atom.name())
                            || metadata.is_some_and(|metadata| {
                                metadata.has_bond(
                                    &endpoint_for_entity(entity),
                                    &BondEndpoint::from_parts(
                                        entity.chain().id(),
                                        entity.residue().serial_number(),
                                        entity.residue().insertion_code().unwrap_or(""),
                                        atom.name(),
                                    ),
                                )
                            }))
                })
        })
        .count()
}

fn is_hydrogen_for_carbon(donor_name: &str, hydrogen_name: &str) -> bool {
    donor_name
        .strip_prefix('C')
        .zip(hydrogen_name.strip_prefix('H'))
        .is_some_and(|(donor_suffix, hydrogen_suffix)| hydrogen_suffix.starts_with(donor_suffix))
}

/// Determine if the two entities are a valid weak hydrogen bond donor-acceptor pair
fn is_weak_donor_acceptor_pair<'a>(
    entity1: &'a AtomConformerResidueChainModel<'a>,
    entity2: &'a AtomConformerResidueChainModel<'a>,
) -> Option<(
    &'a AtomConformerResidueChainModel<'a>,
    &'a AtomConformerResidueChainModel<'a>,
)> {
    let e1_conformer = entity1.conformer().name();
    let e2_conformer = entity2.conformer().name();
    let e1_atom = entity1.atom().name();
    let e2_atom = entity2.atom().name();

    if is_weak_hydrogen_donor(entity1.atom()) && is_hydrogen_acceptor(e2_conformer, e2_atom) {
        Some((entity1, entity2))
    } else if is_weak_hydrogen_donor(entity2.atom()) && is_hydrogen_acceptor(e1_conformer, e1_atom)
    {
        Some((entity2, entity1))
    } else {
        None
    }
}

/// Determine if the atom is a weak hydrogen donor
fn is_weak_hydrogen_donor(atom: &Atom) -> bool {
    // All the non-carbonyl carbon atoms
    (atom.element() == Some(&Element::C)) && atom.name() != "C"
}
