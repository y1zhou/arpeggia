use super::ionic::is_pos_ionizable;
use super::structs::Interaction;
use crate::residues::Ring;

use nalgebra as na;
use pdbtbx::*;

const CATION_PI_ANGLE_THRESHOLD: f64 = 30.0;
const CATION_PI_DIST_THRESHOLD: f64 = 4.5;
const PI_PI_DIST_THRESHOLD: f64 = 6.0;

/// Identify cation-pi interactions.
pub fn find_cation_pi(ring: &Ring, entity: &AtomConformerResidueChainModel) -> Option<Interaction> {
    if is_pos_ionizable(entity.residue().name().unwrap(), entity.atom().name()) {
        let atom_coord = entity.atom().pos();
        let atom_point = na::Vector3::new(atom_coord.0, atom_coord.1, atom_coord.2);
        let dist = point_ring_dist(ring, &atom_coord);
        let theta = point_ring_angle(ring, &atom_point);

        if (theta <= CATION_PI_ANGLE_THRESHOLD) & (dist <= CATION_PI_DIST_THRESHOLD) {
            Some(Interaction::CationPi)
        } else {
            None
        }
    } else {
        None
    }
}

/// Identify pi-pi interactions using the classification in [CREDO](https://doi.org/10.1111/j.1747-0285.2008.00762.x).
pub fn find_pi_pi(ring1: &Ring, ring2: &Ring) -> Option<Interaction> {
    let angle_vec = ring1.center - ring2.center;
    let dist = (angle_vec).norm();
    if dist > PI_PI_DIST_THRESHOLD {
        None
    } else {
        let dihedral = ring_ring_angle(ring1, ring2);
        let theta = point_ring_angle(ring1, &angle_vec);

        match (dihedral, theta) {
            (x, _) if (x > 30.0) => Some(Interaction::PiTStacking),
            (_, y) if (y > 20.0) => Some(Interaction::PiDisplacedStacking),
            _ => Some(Interaction::PiSandwichStacking),
        }
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
