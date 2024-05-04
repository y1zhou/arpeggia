use super::structs::Interaction;

use pdbtbx::*;

/// Search for steric clashes, disulfides, and Van de Waals contacts.
///
/// ## Details
/// If the distance between two atoms is less than the sum of their covalent radii minus `vdw_comp_factor`,
/// they are considered to be in a steric clash.
/// If the distance is between the covalent radii sum and that plus the `vdw_comp_factor`,
/// they are considered to be covalently bonded.
/// A special case would be disulfides, where we also consider the Cb-S-S-Cb dihedral angle
/// and determine if it is around 90° ([10.1039/c8sc01423j](https://doi.org/10.1039/c8sc01423j)).
/// If the distance between two atoms is less than the sum of their Van de Waals radii
/// plus the `vdw_comp_factor`,
/// they are considered to be in a Van de Waals contact.
///
/// TODO: need bonding information for correct covalent bond identification
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
        d if d < sum_cov_radii - vdw_comp_factor => Some(Interaction::StericClash),
        d if d < sum_cov_radii + vdw_comp_factor => match is_disulfide(entity1, entity2) {
            true => Some(Interaction::Disulfide),
            false => Some(Interaction::CovalentBond),
        },
        d if d < sum_vdw_radii + vdw_comp_factor => Some(Interaction::VanDerWaalsContact),
        _ => None,
    }
}

fn is_disulfide(
    entity1: &AtomConformerResidueChainModel,
    entity2: &AtomConformerResidueChainModel,
) -> bool {
    if (entity1.residue().name().unwrap() == "CYS")
        & (entity2.residue().name().unwrap() == "CYS")
        & (entity1.atom().name() == "SG")
        & (entity2.atom().name() == "SG")
    {
        let cb1 = entity1
            .residue()
            .atoms()
            .find(|atom| atom.name() == "CB")
            .unwrap();
        let s1 = entity1
            .residue()
            .atoms()
            .find(|atom| atom.name() == "SG")
            .unwrap();
        let cb2 = entity2
            .residue()
            .atoms()
            .find(|atom| atom.name() == "CB")
            .unwrap();
        let s2 = entity2
            .residue()
            .atoms()
            .find(|atom| atom.name() == "SG")
            .unwrap();
        let cssb_dihedral = cb1.dihedral(s1, s2, cb2).abs();
        (60.0..=120.0).contains(&cssb_dihedral)
    } else {
        false
    }
}
