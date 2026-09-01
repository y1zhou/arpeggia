use super::structs::{Interaction, ProtonationMode};

use pdbtbx::*;

const IONIC_BOND_DIST: f64 = 4.0;

/// Search for ionic bonds.
///
/// Check if the distance between a positive and a negative ionizable atom
/// is within [`IONIC_BOND_DIST`].
/// Search for ionic bonds under an explicit histidine-protonation policy.
pub(super) fn find_ionic_bond_with_protonation(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    mode: ProtonationMode,
    ph: f64,
) -> Option<Interaction> {
    if let Some((pos, neg, charge)) = is_ionic_bond_pair(entity1, entity2, mode, ph)
        && pos.atom().distance(neg.atom()) <= IONIC_BOND_DIST
    {
        return Some(match charge {
            PositiveCharge::Definite => Interaction::IonicBond,
            PositiveCharge::Potential => Interaction::PotentialIonicBond,
        });
    }
    None
}

/// Search for like charges that repel each other.
/// Search for like-charge contacts under a histidine-protonation policy.
pub(super) fn find_ionic_repulsion_with_protonation(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    mode: ProtonationMode,
    ph: f64,
) -> Option<Interaction> {
    if let Some((hier1, hier2, potential)) = is_same_charge_pair(entity1, entity2, mode, ph)
        && hier1.atom().distance(hier2.atom()) <= IONIC_BOND_DIST
    {
        return Some(if potential {
            Interaction::PotentialIonicRepulsion
        } else {
            Interaction::IonicRepulsion
        });
    }
    None
}

/// Check if the two entities are a pair of ionizables.
fn is_ionic_bond_pair<'a>(
    entity1: &'a AtomConformerResidueChainModel<'a>,
    entity2: &'a AtomConformerResidueChainModel<'a>,
    mode: ProtonationMode,
    ph: f64,
) -> Option<(
    &'a AtomConformerResidueChainModel<'a>,
    &'a AtomConformerResidueChainModel<'a>,
    PositiveCharge,
)> {
    let e1_conformer = entity1.conformer().name();
    let e2_conformer = entity2.conformer().name();
    let e1_atom = entity1.atom().name();
    let e2_atom = entity2.atom().name();

    if let Some(charge) = positive_charge(entity1, mode, ph)
        && is_neg_ionizable(e2_conformer, e2_atom)
    {
        Some((entity1, entity2, charge))
    } else if let Some(charge) = positive_charge(entity2, mode, ph)
        && is_neg_ionizable(e1_conformer, e1_atom)
    {
        Some((entity2, entity1, charge))
    } else {
        None
    }
}

fn is_same_charge_pair<'a>(
    entity1: &'a AtomConformerResidueChainModel<'a>,
    entity2: &'a AtomConformerResidueChainModel<'a>,
    mode: ProtonationMode,
    ph: f64,
) -> Option<(
    &'a AtomConformerResidueChainModel<'a>,
    &'a AtomConformerResidueChainModel<'a>,
    bool,
)> {
    let e1_conformer = entity1.conformer().name();
    let e2_conformer = entity2.conformer().name();
    let e1_atom = entity1.atom().name();
    let e2_atom = entity2.atom().name();

    let charge1 = positive_charge(entity1, mode, ph);
    let charge2 = positive_charge(entity2, mode, ph);
    let both_pos_charged = charge1.is_some() && charge2.is_some();
    let both_neg_charged =
        is_neg_ionizable(e1_conformer, e1_atom) && is_neg_ionizable(e2_conformer, e2_atom);

    if both_pos_charged | both_neg_charged {
        Some((
            entity1,
            entity2,
            charge1 == Some(PositiveCharge::Potential)
                || charge2 == Some(PositiveCharge::Potential),
        ))
    } else {
        None
    }
}

/// Check if the entity contains ionizable groups that are positively charged at pH 7.0.
fn is_pos_ionizable(res_name: &str, atom_name: &str) -> bool {
    matches!(
        (res_name, atom_name),
        ("ARG", "NE" | "CZ" | "NH1" | "NH2")
            | (
                "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP",
                "CG" | "ND1" | "CE1" | "NE2" | "CD2"
            )
            | ("LYS", "NZ")
    )
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum PositiveCharge {
    Definite,
    Potential,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum HistidinePreparationIssue {
    Unresolved,
    Inconsistent,
}

pub(super) fn positive_charge(
    entity: &AtomConformerResidueChainModel,
    mode: ProtonationMode,
    ph: f64,
) -> Option<PositiveCharge> {
    let residue_name = entity.conformer().name();
    let atom_name = entity.atom().name();
    if !is_histidine(residue_name) {
        return is_pos_ionizable(residue_name, atom_name).then_some(PositiveCharge::Definite);
    }
    if !is_pos_ionizable(residue_name, atom_name) {
        return None;
    }
    let explicit = explicit_histidine_charge(entity.residue());
    if explicit == Some(true) {
        return Some(PositiveCharge::Definite);
    }
    if mode == ProtonationMode::AllCharged {
        return Some(PositiveCharge::Potential);
    }
    match explicit {
        Some(false) => None,
        None => match mode {
            ProtonationMode::AllCharged => unreachable!(),
            // [WARNING] An intrinsic-pKa threshold is a population prior, not a
            // per-residue structure-based protonation prediction.
            ProtonationMode::Heuristic => (ph <= 6.0).then_some(PositiveCharge::Potential),
            ProtonationMode::ExplicitOnly => None,
        },
        Some(true) => unreachable!(),
    }
}

pub(crate) fn is_histidine(name: &str) -> bool {
    matches!(name, "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP")
}

fn explicit_histidine_charge(residue: &Residue) -> Option<bool> {
    if residue.atoms().map(Atom::charge).sum::<isize>() > 0 {
        return Some(true);
    }
    match residue.name().unwrap_or("") {
        "HSP" => return Some(true),
        "HID" | "HIE" | "HSD" | "HSE" => return Some(false),
        "HIP"
            if residue
                .atoms()
                .all(|atom| atom.element() != Some(&Element::P)) =>
        {
            return Some(true);
        }
        _ => {}
    }
    let hd1 = residue.atoms().any(|atom| atom.name() == "HD1");
    let he2 = residue.atoms().any(|atom| atom.name() == "HE2");
    match (hd1, he2) {
        (true, true) => Some(true),
        (true, false) | (false, true) => Some(false),
        (false, false) => None,
    }
}

pub(crate) fn histidine_preparation_issue(residue: &Residue) -> Option<HistidinePreparationIssue> {
    let name = residue.name().unwrap_or("");
    if !is_histidine(name) {
        return None;
    }
    let hd1 = residue.atoms().any(|atom| atom.name() == "HD1");
    let he2 = residue.atoms().any(|atom| atom.name() == "HE2");
    let inconsistent = match name {
        "HID" | "HSD" => !hd1 || he2,
        "HIE" | "HSE" => hd1 || !he2,
        "HIP" | "HSP" => {
            !hd1 || !he2
                || residue
                    .atoms()
                    .any(|atom| atom.element() == Some(&Element::P))
        }
        _ => false,
    };
    if inconsistent {
        Some(HistidinePreparationIssue::Inconsistent)
    } else if name == "HIS" && !hd1 && !he2 {
        Some(HistidinePreparationIssue::Unresolved)
    } else {
        None
    }
}

/// Check if the entity contains ionizable groups that are negatively charged at pH 7.0.
fn is_neg_ionizable(res_name: &str, atom_name: &str) -> bool {
    matches!(
        (res_name, atom_name),
        ("ASP", "OD1" | "OD2") | ("GLU", "OE1" | "OE2")
    )
}
