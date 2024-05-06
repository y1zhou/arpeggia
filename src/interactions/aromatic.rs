use super::ionic::is_pos_ionizable;
use super::structs::Interaction;
use crate::residues::Ring;

use nalgebra as na;
use pdbtbx::*;

pub fn find_cation_pi(ring: &Ring, entity: &AtomConformerResidueChainModel) -> Option<Interaction> {
    if is_pos_ionizable(entity.residue().name().unwrap(), entity.atom().name()) {
        let atom_coord = entity.atom().pos();
        let atom_point = na::Vector3::new(atom_coord.0, atom_coord.1, atom_coord.2);
        let dist = point_ring_dist(ring, &atom_coord);
        let theta = point_ring_angle(ring, &atom_point);

        if (theta <= 30.0) & (dist <= 4.5) {
            Some(Interaction::CationPi)
        } else {
            None
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
    let mut rad = (ring.normal.dot(&v) / (ring.normal.norm() * point.norm())).acos();

    // Convert to degrees
    if rad > std::f64::consts::FRAC_PI_2 {
        rad -= std::f64::consts::PI;
    }
    rad.to_degrees()
}

/// Calculate the angle between two rings.
fn ring_ring_angle(ring1: &Ring, ring2: &Ring) -> f64 {
    point_ring_angle(ring1, &ring2.normal)
}
