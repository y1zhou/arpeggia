use super::structs::Interaction;
use crate::metadata::ResolvedBonds;

use pdbtbx::*;

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
pub(super) fn find_hydrogen_bond(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
    bonds: &ResolvedBonds,
) -> Option<Interaction> {
    if let Some((donor, acceptor)) = is_donor_acceptor_pair(entity1, entity2) {
        let da_dist = donor.atom().distance(acceptor.atom());
        if da_dist <= HYDROGEN_BOND_DIST {
            let donor_has_valid_h = donor
                .residue()
                .atoms()
                .filter(|atom| {
                    atom.element() == Some(&Element::H)
                        && (is_hydrogen_for_donor(
                            donor.conformer().name(),
                            donor.atom().name(),
                            atom.name(),
                        ) || bonds.contains_entity_atom(donor, atom))
                })
                .any(|hydrogen| {
                    let acceptor_vdw = acceptor
                        .atom()
                        .element()
                        .and_then(|element| element.atomic_radius().van_der_waals);
                    let h_vdw = Element::H.atomic_radius().van_der_waals;
                    acceptor_vdw
                        .zip(h_vdw)
                        .is_some_and(|(acceptor_vdw, h_vdw)| {
                            hydrogen.distance(acceptor.atom())
                                <= h_vdw + acceptor_vdw + vdw_comp_factor
                                && donor.atom().angle(hydrogen, acceptor.atom()) >= 90.0
                        })
                });
            if donor_has_valid_h {
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
pub(super) fn find_weak_hydrogen_bond(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
    bonds: &ResolvedBonds,
) -> Option<Interaction> {
    if let Some((donor, acceptor)) = is_weak_donor_acceptor_pair(entity1, entity2) {
        let da_dist = donor.atom().distance(acceptor.atom());
        if da_dist <= HYDROGEN_BOND_DIST {
            let donor_has_valid_h = donor
                .residue()
                .atoms()
                .filter(|atom| {
                    atom.element() == Some(&Element::H)
                        && (is_hydrogen_for_carbon(donor.atom().name(), atom.name())
                            || bonds.contains_entity_atom(donor, atom))
                })
                .any(|hydrogen| {
                    let acceptor_vdw = acceptor
                        .atom()
                        .element()
                        .and_then(|element| element.atomic_radius().van_der_waals);
                    let h_vdw = Element::H.atomic_radius().van_der_waals;
                    acceptor_vdw
                        .zip(h_vdw)
                        .is_some_and(|(acceptor_vdw, h_vdw)| {
                            hydrogen.distance(acceptor.atom())
                                <= h_vdw + acceptor_vdw + vdw_comp_factor
                                && donor.atom().angle(hydrogen, acceptor.atom()) >= 130.0
                        })
                });
            if donor_has_valid_h {
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

    if is_hydrogen_donor(e1_conformer, e1_atom, entity1.residue())
        && is_hydrogen_acceptor(e2_conformer, e2_atom, entity2.residue())
    {
        Some((entity1, entity2))
    } else if is_hydrogen_donor(e2_conformer, e2_atom, entity2.residue())
        && is_hydrogen_acceptor(e1_conformer, e1_atom, entity1.residue())
    {
        Some((entity2, entity1))
    } else {
        None
    }
}

/// Determine if the atom in the residue is a hydrogen acceptor
fn is_hydrogen_acceptor(res_name: &str, atom_name: &str, residue: &Residue) -> bool {
    // all the carbonyl oxygens in the main chain and on the terminals
    if matches!(atom_name, "O" | "OXT") && res_name != "HOH" {
        return true;
    }
    if crate::contacts::ionic::is_histidine(res_name) && matches!(atom_name, "ND1" | "NE2") {
        return match histidine_tautomer(res_name, residue) {
            HistidineTautomer::Delta => atom_name == "NE2",
            HistidineTautomer::Epsilon => atom_name == "ND1",
            HistidineTautomer::Both => false,
            HistidineTautomer::Unresolved => true,
        };
    }
    matches!(
        (res_name, atom_name),
        ("ASN", "OD1")
            // | ("ASN", "ND2")
            | ("ASP", "OD1" | "OD2")
            | ("GLN", "OE1")
            | ("GLU", "OE1" | "OE2")
            | ("SER", "OG")
            | ("THR", "OG1")
            | ("TYR", "OH")
            | ("MET", "SD") // 10.1021/jz300207k and 10.1002/prot.22327
            | ("CYS", "SG") // 10.1002/prot.22327
    )
}

/// Determine if the atom in the residue is a hydrogen donor
fn is_hydrogen_donor(res_name: &str, atom_name: &str, residue: &Residue) -> bool {
    // Backbone amide nitrogens donate; proline does so only when the input
    // explicitly supplies a terminal N-bound hydrogen.
    if atom_name == "N" {
        return res_name != "PRO"
            || residue.atoms().any(|atom| {
                atom.element() == Some(&Element::H)
                    && is_hydrogen_for_donor(res_name, atom_name, atom.name())
            });
    }
    if crate::contacts::ionic::is_histidine(res_name) && matches!(atom_name, "ND1" | "NE2") {
        return match histidine_tautomer(res_name, residue) {
            HistidineTautomer::Delta => atom_name == "ND1",
            HistidineTautomer::Epsilon => atom_name == "NE2",
            HistidineTautomer::Both | HistidineTautomer::Unresolved => true,
        };
    }
    matches!(
        (res_name, atom_name),
        ("ARG", "NE" | "NH1" | "NH2")
            | ("ASN", "ND2")
            | ("GLN", "NE2")
            | ("LYS", "NZ")
            | ("SER", "OG")
            | ("THR", "OG1")
            | ("TRP", "NE1")
            | ("TYR", "OH")
            | ("CYS", "SG") // 10.1002/prot.22327
    )
}

#[derive(Clone, Copy)]
enum HistidineTautomer {
    Delta,
    Epsilon,
    Both,
    Unresolved,
}

fn histidine_tautomer(res_name: &str, residue: &Residue) -> HistidineTautomer {
    match res_name {
        "HID" | "HSD" => HistidineTautomer::Delta,
        "HIE" | "HSE" => HistidineTautomer::Epsilon,
        "HIP" | "HSP" => HistidineTautomer::Both,
        _ => {
            let hd1 = residue.atoms().any(|atom| atom.name() == "HD1");
            let he2 = residue.atoms().any(|atom| atom.name() == "HE2");
            match (hd1, he2) {
                (true, false) => HistidineTautomer::Delta,
                (false, true) => HistidineTautomer::Epsilon,
                (true, true) => HistidineTautomer::Both,
                (false, false) => HistidineTautomer::Unresolved,
            }
        }
    }
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

pub(super) fn count_donors_without_explicit_hydrogen(
    pdb: &PDB,
    bonds: &ResolvedBonds,
    ligand: &std::collections::HashSet<String>,
    receptor: &std::collections::HashSet<String>,
) -> usize {
    pdb.atoms_with_hierarchy()
        .filter(|entity| {
            ligand.contains(entity.chain().id()) || receptor.contains(entity.chain().id())
        })
        .filter(|entity| {
            let residue_name = entity.conformer().name();
            let donor_name = entity.atom().name();
            is_hydrogen_donor(residue_name, donor_name, entity.residue())
                && !entity.residue().atoms().any(|atom| {
                    atom.element() == Some(&Element::H)
                        && (is_hydrogen_for_donor(residue_name, donor_name, atom.name())
                            || bonds.contains_entity_atom(entity, atom))
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

    if is_weak_hydrogen_donor(e1_conformer, entity1.atom())
        && is_hydrogen_acceptor(e2_conformer, e2_atom, entity2.residue())
    {
        Some((entity1, entity2))
    } else if is_weak_hydrogen_donor(e2_conformer, entity2.atom())
        && is_hydrogen_acceptor(e1_conformer, e1_atom, entity1.residue())
    {
        Some((entity2, entity1))
    } else {
        None
    }
}

/// Determine if the atom is a weak hydrogen donor
fn is_weak_hydrogen_donor(residue: &str, atom: &Atom) -> bool {
    atom.element() == Some(&Element::C)
        && matches!(
            (residue, atom.name()),
            (
                "ALA"
                    | "ARG"
                    | "ASN"
                    | "ASP"
                    | "CYS"
                    | "GLN"
                    | "GLU"
                    | "GLY"
                    | "HIS"
                    | "HID"
                    | "HIE"
                    | "HIP"
                    | "HSD"
                    | "HSE"
                    | "HSP"
                    | "ILE"
                    | "LEU"
                    | "LYS"
                    | "MET"
                    | "MSE"
                    | "PHE"
                    | "PRO"
                    | "SER"
                    | "THR"
                    | "TRP"
                    | "TYR"
                    | "VAL",
                "CA"
            ) | ("ALA", "CB")
                | ("ARG", "CB" | "CG" | "CD")
                | ("ASN" | "ASP" | "CYS" | "SER", "CB")
                | ("GLN" | "GLU", "CB" | "CG")
                | (
                    "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP",
                    "CB" | "CD2" | "CE1"
                )
                | ("ILE", "CB" | "CG1" | "CG2" | "CD1")
                | ("LEU", "CB" | "CG" | "CD1" | "CD2")
                | ("LYS", "CB" | "CG" | "CD" | "CE")
                | ("MET" | "MSE", "CB" | "CG" | "CE")
                | ("PHE", "CB" | "CD1" | "CD2" | "CE1" | "CE2" | "CZ")
                | ("PRO", "CB" | "CG" | "CD")
                | ("THR", "CB" | "CG2")
                | ("TRP", "CB" | "CD1" | "CE3" | "CZ2" | "CZ3" | "CH2")
                | ("TYR", "CB" | "CD1" | "CD2" | "CE1" | "CE2")
                | ("VAL", "CB" | "CG1" | "CG2")
        )
}
