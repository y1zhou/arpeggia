use super::structs::Interaction;

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
/// If so, it collects all the hydrogen atoms of the donor residue and checks if any of them satisfy:
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
/// TODO: add special case for water-mediated hydrogen bonds
pub fn find_hydrogen_bond(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
) -> Option<Interaction> {
    if let Some((donor, acceptor)) = is_donor_acceptor_pair(entity1, entity2) {
        let da_dist = donor.atom().distance(acceptor.atom());
        if da_dist <= HYDROGEN_BOND_DIST {
            let donor_h: Vec<&Atom> = donor
                .residue()
                .par_atoms()
                .filter(|atom| atom.element().unwrap() == &Element::H)
                .collect();

            // Hydrogen bonds are stricter as they have angle restrictions
            let acceptor_vdw: f64 = acceptor
                .atom()
                .element()
                .unwrap()
                .atomic_radius()
                .van_der_waals
                .unwrap();
            let h_vdw: f64 = Element::H.atomic_radius().van_der_waals.unwrap();
            if donor_h.par_iter().any(|h| {
                (h.distance(acceptor.atom()) <= h_vdw + acceptor_vdw + vdw_comp_factor)
                    & (donor.atom().angle(h, acceptor.atom()) >= 90.0)
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
pub fn find_weak_hydrogen_bond(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
) -> Option<Interaction> {
    if let Some((donor, acceptor)) = is_weak_donor_acceptor_pair(entity1, entity2) {
        let da_dist = donor.atom().distance(acceptor.atom());
        if da_dist <= HYDROGEN_BOND_DIST {
            let donor_h: Vec<&Atom> = donor
                .residue()
                .par_atoms()
                .filter(|atom| atom.element().unwrap() == &Element::H)
                .collect();

            // Hydrogen bonds are stricter as they have angle restrictions
            let acceptor_vdw: f64 = acceptor
                .atom()
                .element()
                .unwrap()
                .atomic_radius()
                .van_der_waals
                .unwrap();
            let h_vdw: f64 = Element::H.atomic_radius().van_der_waals.unwrap();
            if donor_h.par_iter().any(|h| {
                (h.distance(acceptor.atom()) <= h_vdw + acceptor_vdw + vdw_comp_factor)
                    & (donor.atom().angle(h, acceptor.atom()) >= 130.0)
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

    if is_hydrogen_donor(e1_conformer, e1_atom) & is_hydrogen_acceptor(e2_conformer, e2_atom) {
        Some((entity1, entity2))
    } else if is_hydrogen_donor(e2_conformer, e2_atom) & is_hydrogen_acceptor(e1_conformer, e1_atom)
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
            | ("ASP", "OD1")
            | ("ASP", "OD2")
            | ("GLN", "OE1")
            | ("GLU", "OE1")
            | ("GLU", "OE2")
            | ("HIS", "ND1")
            | ("HIS", "NE2")
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
    if atom_name == "N" {
        return true;
    }
    matches!(
        (res_name, atom_name),
        ("ARG", "NE")
            | ("ARG", "NH1")
            | ("ARG", "NH2")
            | ("ASN", "ND2")
            | ("GLN", "NE2")
            | ("HIS", "ND1")
            | ("HIS", "NE2")
            | ("LYS", "NZ")
            | ("SER", "OG")
            | ("THR", "OG1")
            | ("TRP", "NE1")
            | ("TYR", "OH")
            | ("CYS", "SG") // 10.1002/prot.22327
    )
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

    if is_weak_hydrogen_donor(entity1.atom()) & is_hydrogen_acceptor(e2_conformer, e2_atom) {
        Some((entity1, entity2))
    } else if is_weak_hydrogen_donor(entity2.atom()) & is_hydrogen_acceptor(e1_conformer, e1_atom) {
        Some((entity2, entity1))
    } else {
        None
    }
}

/// Determine if the atom is a weak hydrogen donor
fn is_weak_hydrogen_donor(atom: &Atom) -> bool {
    // All the non-carbonyl carbon atoms
    (atom.element().unwrap() == &Element::C) && atom.name() != "C"
}
