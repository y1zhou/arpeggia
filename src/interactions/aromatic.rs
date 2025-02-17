use super::ionic::is_pos_ionizable;
use super::structs::Interaction;
use crate::residues::Ring;

use nalgebra as na;
use pdbtbx::*;

const CATION_PI_ANGLE_THRESHOLD: f64 = 30.0;
const CATION_PI_DIST_THRESHOLD: f64 = 4.5;
const PI_PI_DIST_THRESHOLD: f64 = 6.0;
const PI_T_DIST_THREHOLD: f64 = 5.0;

/// Identify cation-pi interactions.
pub fn find_cation_pi(ring: &Ring, entity: &AtomConformerResidueChainModel) -> Option<Interaction> {
    if is_pos_ionizable(entity.residue().name().unwrap(), entity.atom().name()) {
        let atom_coord = entity.atom().pos();
        let atom_point = na::Vector3::new(atom_coord.0, atom_coord.1, atom_coord.2);
        let dist = point_ring_dist(ring, &atom_coord);
        let theta = point_ring_angle(ring, &atom_point);

        if (theta <= CATION_PI_ANGLE_THRESHOLD) & (dist <= CATION_PI_DIST_THRESHOLD) {
            return Some(Interaction::CationPi);
        }
    }
    None
}

/// Identify pi-pi interactions using the classification by [Chakrabarti and Bhattacharyya (2007)](https://doi.org/10.1016/j.pbiomolbio.2007.03.016), Fig. 11.
/// For T-shaped Pi-stacking, the distance threshold is set to 5.0 Å according to [getcontacts](https://getcontacts.github.io/interactions.html).
pub fn find_pi_pi(ring1: &Ring, ring2: &Ring) -> Option<Interaction> {
    let angle_vec = ring1.center - ring2.center;
    let dist = (angle_vec).norm();
    if dist <= PI_PI_DIST_THRESHOLD {
        let theta = point_ring_angle(ring1, &ring2.center);
        let dihedral = ring_ring_angle(ring1, ring2);

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

/// Calculate the distance from the point to the ring center.
pub fn point_ring_dist(ring: &Ring, point: &(f64, f64, f64)) -> f64 {
    let atom_point = na::Vector3::new(point.0, point.1, point.2);
    (atom_point - ring.center).norm()
}

/// Calculate the angle between the ring normal and the vector pointing from
/// the ring center to the point.
fn point_ring_angle(ring: &Ring, point: &na::Vector3<f64>) -> f64 {
    let v = point - ring.center;
    let mut rad = (ring.normal.dot(&v) / (ring.normal.norm() * v.norm())).acos();

    // Convert to degrees
    if rad > std::f64::consts::FRAC_PI_2 {
        rad = std::f64::consts::PI - rad;
    }
    rad.to_degrees()
}

/// Calculate the angle between two rings.
fn ring_ring_angle(ring1: &Ring, ring2: &Ring) -> f64 {
    let mut rad =
        (ring1.normal.dot(&ring2.normal) / (ring1.normal.norm() * ring2.normal.norm())).acos();

    // Convert to degrees
    if rad > std::f64::consts::FRAC_PI_2 {
        rad = std::f64::consts::PI - rad;
    }
    rad.to_degrees()
}
