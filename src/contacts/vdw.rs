use super::structs::Interaction;
use crate::metadata::ExplicitBondKind;

use pdbtbx::*;

/// Search for steric clashes, potential covalency, and Van de Waals contacts.
///
/// ## Details
/// If the distance between two atoms is less than the sum of their covalent radii minus `vdw_comp_factor`,
/// they are considered to be in a steric clash.
/// If the distance is between the covalent radii sum and that plus the `vdw_comp_factor`,
/// they are considered potentially covalent until input topology confirms a bond.
/// Cysteine SG pairs in this band are classified as potential disulfides when
/// their CB-SG-SG-CB dihedral is between 60 and 120 degrees, following the
/// original Arpeggia rule ([DOI: 10.1039/c8sc01423j](https://doi.org/10.1039/c8sc01423j)).
/// Resolved disulfide declarations take precedence over other explicit
/// connectivity, which in turn takes precedence over geometry inference.
/// Nonbonded atoms inside their unmodified Van de Waals envelope clash; atoms only
/// inside the compensated outer envelope are Van de Waals contacts.
pub(super) fn find_vdw_contact(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
    explicit_bond: Option<ExplicitBondKind>,
) -> Option<Interaction> {
    match explicit_bond {
        Some(ExplicitBondKind::Disulfide) => return Some(Interaction::Disulfide),
        Some(ExplicitBondKind::Covalent) => return Some(Interaction::Covalent),
        None => {}
    }
    let e1_radii = entity1.atom().element()?.atomic_radius();
    let e2_radii = entity2.atom().element()?.atomic_radius();

    let sum_cov_radii = e1_radii.covalent_single + e2_radii.covalent_single;
    let sum_vdw_radii = e1_radii.van_der_waals? + e2_radii.van_der_waals?;

    let dist = entity1.atom().distance(entity2.atom());

    match dist {
        d if d < sum_cov_radii - vdw_comp_factor => Some(Interaction::StericClash),
        d if d < sum_cov_radii + vdw_comp_factor && is_disulfide(entity1, entity2) => {
            Some(Interaction::PotentialDisulfide)
        }
        d if d < sum_cov_radii + vdw_comp_factor => Some(Interaction::PotentialCovalent),
        d if d < sum_vdw_radii => Some(Interaction::VanDerWaalsClash),
        d if d < sum_vdw_radii + vdw_comp_factor => Some(Interaction::VanDerWaalsContact),
        _ => None,
    }
}

fn is_disulfide(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
) -> bool {
    if entity1.residue().name() != Some("CYS")
        || entity2.residue().name() != Some("CYS")
        || entity1.atom().name() != "SG"
        || entity2.atom().name() != "SG"
    {
        return false;
    }

    let Some(cb1) = entity1.residue().atoms().find(|atom| atom.name() == "CB") else {
        return false;
    };
    let Some(cb2) = entity2.residue().atoms().find(|atom| atom.name() == "CB") else {
        return false;
    };

    (60.0..=120.0).contains(&cb1.dihedral(entity1.atom(), entity2.atom(), cb2).abs())
}
