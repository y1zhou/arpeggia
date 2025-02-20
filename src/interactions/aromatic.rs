use super::ionic::is_pos_ionizable;
use super::structs::Interaction;
use crate::residues::Plane;

use nalgebra as na;
use pdbtbx::*;

const CATION_PI_ANGLE_THRESHOLD: f64 = 30.0;
const CATION_PI_DIST_THRESHOLD: f64 = 4.5;
const PI_PI_DIST_THRESHOLD: f64 = 6.0;
const PI_T_DIST_THREHOLD: f64 = 5.0;

/// Identify cation-pi interactions.
pub fn find_cation_pi(
    ring: &Plane,
    entity: &AtomConformerResidueChainModel,
) -> Option<Interaction> {
    if is_pos_ionizable(entity.residue().name().unwrap(), entity.atom().name()) {
        let atom_coord = entity.atom().pos();
        let atom_point = na::Vector3::new(atom_coord.0, atom_coord.1, atom_coord.2);
        let dist = ring.point_vec_dist(&atom_point);
        let theta = ring.point_vec_angle(&atom_point);

        if (theta <= CATION_PI_ANGLE_THRESHOLD) & (dist <= CATION_PI_DIST_THRESHOLD) {
            return Some(Interaction::CationPi);
        }
    }
    None
}

/// Identify pi-pi interactions using the classification by [Chakrabarti and Bhattacharyya (2007)](https://doi.org/10.1016/j.pbiomolbio.2007.03.016), Fig. 11.
/// For T-shaped Pi-stacking, the distance threshold is set to 5.0 Å according to [getcontacts](https://getcontacts.github.io/interactions.html).
pub fn find_pi_pi(ring1: &Plane, ring2: &Plane) -> Option<Interaction> {
    let angle_vec = ring1.center - ring2.center;
    let dist = (angle_vec).norm();
    if dist <= PI_PI_DIST_THRESHOLD {
        let theta = ring1.point_vec_angle(&ring2.center);
        let dihedral = ring1.dihedral(ring2);

        match dihedral {
            d if d <= 30.0 => match theta {
                t if t <= 30.0 => Some(Interaction::PiSandwichStacking), // ff
                t if t <= 60.0 => Some(Interaction::PiDisplacedStacking), // of
                t if t <= 90.0 => Some(Interaction::PiParallelInPlaneStacking), // ee
                _ => None,
            },
            // ft, ot, and et
            d if d <= 60.0 => Some(Interaction::PiTiltedStacking),
            d if d <= 90.0 => match theta {
                t if (30.0..60.0).contains(&t) => Some(Interaction::PiLStacking), // oe
                _ => match dist <= PI_T_DIST_THREHOLD {
                    true => Some(Interaction::PiTStacking), // fe and ef
                    false => None,
                },
            },
            _ => None,
        }
    } else {
        None
    }
}
