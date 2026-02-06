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
    if let Some((pos, neg)) = is_ionic_bond_pair(entity1, entity2)
        && pos.atom().distance(neg.atom()) <= IONIC_BOND_DIST
    {
        // Polar interactions are more relaxed as only distance is checked
        return Some(Interaction::IonicBond);
    }
    None
}

/// Search for like charges that repel each other.
pub fn find_ionic_repulsion(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
) -> Option<Interaction> {
    if let Some((hier1, hier2)) = is_same_charge_pair(entity1, entity2)
        && hier1.atom().distance(hier2.atom()) <= IONIC_BOND_DIST
    {
        return Some(Interaction::IonicRepulsion);
    }
    None
}

/// Check if the two entities are a pair of ionizables.
fn is_ionic_bond_pair<'a>(
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

    if is_pos_ionizable(e1_conformer, e1_atom) && is_neg_ionizable(e2_conformer, e2_atom) {
        Some((entity1, entity2))
    } else if is_pos_ionizable(e2_conformer, e2_atom) & is_neg_ionizable(e1_conformer, e1_atom) {
        Some((entity2, entity1))
    } else {
        None
    }
}

fn is_same_charge_pair<'a>(
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

    let both_pos_charged =
        is_pos_ionizable(e1_conformer, e1_atom) && is_pos_ionizable(e2_conformer, e2_atom);
    let both_neg_charged =
        is_neg_ionizable(e1_conformer, e1_atom) && is_neg_ionizable(e2_conformer, e2_atom);

    if both_pos_charged | both_neg_charged {
        Some((entity1, entity2))
    } else {
        None
    }
}

/// Check if the entity contains ionizable groups that are positively charged at pH 7.0.
pub fn is_pos_ionizable(res_name: &str, atom_name: &str) -> bool {
    matches!(
        (res_name, atom_name),
        ("ARG", "NE" | "CZ" | "NH1" | "NH2")
            | ("HIS", "CG" | "ND1" | "CE1" | "NE2" | "CD2")
            | ("LYS", "NZ")
    )
}

/// Check if the entity contains ionizable groups that are negatively charged at pH 7.0.
fn is_neg_ionizable(res_name: &str, atom_name: &str) -> bool {
    matches!(
        (res_name, atom_name),
        ("ASP", "OD1" | "OD2") | ("GLU", "OE1" | "OE2")
    )
}
