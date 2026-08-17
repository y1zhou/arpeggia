use super::structs::Interaction;

use pdbtbx::*;

/// Search for steric clashes, potential covalency, and Van de Waals contacts.
///
/// ## Details
/// If the distance between two atoms is less than the sum of their covalent radii minus `vdw_comp_factor`,
/// they are considered to be in a steric clash.
/// If the distance is between the covalent radii sum and that plus the `vdw_comp_factor`,
/// they are considered potentially covalent until input topology confirms a bond.
/// Nonbonded atoms inside their unmodified Van de Waals envelope clash; atoms only
/// inside the compensated outer envelope are Van de Waals contacts.
pub fn find_vdw_contact(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
    explicitly_bonded: bool,
) -> Option<Interaction> {
    if explicitly_bonded {
        return Some(Interaction::Covalent);
    }
    let e1_radii = entity1.atom().element()?.atomic_radius();
    let e2_radii = entity2.atom().element()?.atomic_radius();

    let sum_cov_radii = e1_radii.covalent_single + e2_radii.covalent_single;
    let sum_vdw_radii = e1_radii.van_der_waals? + e2_radii.van_der_waals?;

    let dist = entity1.atom().distance(entity2.atom());

    match dist {
        d if d < sum_cov_radii - vdw_comp_factor => Some(Interaction::StericClash),
        d if d < sum_cov_radii + vdw_comp_factor => Some(Interaction::PotentialCovalent),
        d if d < sum_vdw_radii => Some(Interaction::VanDerWaalsClash),
        d if d < sum_vdw_radii + vdw_comp_factor => Some(Interaction::VanDerWaalsContact),
        _ => None,
    }
}
