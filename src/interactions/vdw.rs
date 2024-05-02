use super::structs::Interaction;

use pdbtbx::*;

pub fn find_vdw_contact(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
    vdw_comp_factor: f64,
) -> Option<Interaction> {
    let e1_radii = entity1.atom().element().unwrap().atomic_radius();
    let e2_radii = entity2.atom().element().unwrap().atomic_radius();

    let sum_cov_radii = e1_radii.covalent_single + e2_radii.covalent_single;
    let sum_vdw_radii = e1_radii.van_der_waals.unwrap() + e2_radii.van_der_waals.unwrap();

    let dist = entity1.atom().distance(entity2.atom());

    match dist {
        d if d < sum_cov_radii => Some(Interaction::StericClash),
        // TODO: need bonding information for correct covalent bond identification
        // d if d < sum_vdw_radii => Some(Interaction::CovalentBond),
        d if d < sum_vdw_radii + vdw_comp_factor => Some(Interaction::VanDerWaalsContact),
        _ => None,
    }
}
