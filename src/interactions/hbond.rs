use super::structs::Interaction;

use pdbtbx::*;
use rayon::prelude::*;
use std::collections::HashSet;

const HYDROGEN_BOND_POLAR_DIST: f64 = 3.5;

/// Search for hydrogen bonds and polar contacts.
///
/// `vdw_comp_factor` is the compensation factor for VdW radii dependent interaction types.
///
/// See details in:
/// https://pymolwiki.org/index.php/Displaying_Biochemical_Properties#Hydrogen_bonds_and_Polar_Contacts
///
/// TODO: add special case for water-mediated hydrogen bonds
pub fn find_hydrogen_bond(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
) -> Option<Interaction> {
    if let Some((donor, acceptor)) = is_donor_acceptor_pair(entity1, entity2) {
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
            Some(Interaction::HydrogenBond)
        } else if donor.atom().distance(acceptor.atom()) <= HYDROGEN_BOND_POLAR_DIST {
            // Polar interactions are more relaxed as only distance is checked
            Some(Interaction::PolarContact)
        } else {
            None
        }
    } else {
        None
    }
}

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

fn is_weak_hydrogen_donor(atom: &Atom) -> bool {
    // All the non-carbonyl carbon atoms
    (atom.element().unwrap() == &Element::C) && atom.name() != "C"
}
