use super::structs::Interaction;

use pdbtbx::*;

const IONIC_BOND_DIST: f64 = 4.0;

/// Search for ionic bonds.
///
/// Check if the distance between a positive and a negative ionizable atom
/// is within [`IONIC_BOND_DIST`].
pub fn find_ionic_bond(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
) -> Option<Interaction> {
    if let Some((pos, neg)) = is_ionic_pair(entity1, entity2) {
        if pos.atom().distance(neg.atom()) <= IONIC_BOND_DIST {
            // Polar interactions are more relaxed as only distance is checked
            Some(Interaction::IonicBond)
        } else {
            None
        }
    } else {
        None
    }
}

/// Check if the two entities are a pair of ionizables.
fn is_ionic_pair<'a>(
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

    if is_pos_ionizable(e1_conformer, e1_atom) & is_neg_ionizable(e2_conformer, e2_atom) {
        Some((entity1, entity2))
    } else if is_pos_ionizable(e2_conformer, e2_atom) & is_neg_ionizable(e1_conformer, e1_atom) {
        Some((entity2, entity1))
    } else {
        None
    }
}

/// Check if the entity contains ionizable groups that are positively charged at pH 7.0.
fn is_pos_ionizable(res_name: &str, atom_name: &str) -> bool {
    matches!(
        (res_name, atom_name),
        ("ARG", "NE")
            | ("ARG", "CZ")
            | ("ARG", "NH1")
            | ("ARG", "NH2")
            | ("HIS", "CG")
            | ("HIS", "ND1")
            | ("HIS", "CE1")
            | ("HIS", "NE2")
            | ("HIS", "CD2")
            | ("LYS", "NZ")
    )
}

/// Check if the entity contains ionizable groups that are negatively charged at pH 7.0.
fn is_neg_ionizable(res_name: &str, atom_name: &str) -> bool {
    matches!(
        (res_name, atom_name),
        ("ASP", "OD1") | ("ASP", "OD2") | ("GLU", "OE1") | ("GLU", "OE2")
    )
}
