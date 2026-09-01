use super::ProtonationMode;
use super::ionic::{PositiveCharge, positive_charge};
use super::residues::Plane;
use super::structs::Interaction;

use nalgebra as na;
use pdbtbx::*;

const CATION_PI_ANGLE_THRESHOLD: f64 = 30.0;
const CATION_PI_DIST_THRESHOLD: f64 = 4.5;
const PI_PI_DIST_THRESHOLD: f64 = 6.0;
const PI_T_DIST_THREHOLD: f64 = 5.0;

/// Identify cation-pi interactions.
pub(super) fn find_cation_pi(
    ring: &Plane,
    entity: &AtomConformerResidueChainModel,
    protonation: ProtonationMode,
    ph: f64,
) -> Option<Interaction> {
    if let Some(charge) = positive_charge(entity, protonation, ph) {
        let atom_coord = entity.atom().pos();
        let atom_point = na::Vector3::new(atom_coord.0, atom_coord.1, atom_coord.2);
        let dist = ring.point_vec_dist(&atom_point);
        let theta = ring.point_vec_angle(&atom_point);

        if (theta <= CATION_PI_ANGLE_THRESHOLD) && (dist <= CATION_PI_DIST_THRESHOLD) {
            return Some(match charge {
                PositiveCharge::Definite => Interaction::CationPi,
                PositiveCharge::Potential => Interaction::PotentialCationPi,
            });
        }
    }
    None
}

/// Identify pi-pi interactions using the classification by [Chakrabarti and Bhattacharyya (2007)](https://doi.org/10.1016/j.pbiomolbio.2007.03.016), Fig. 11.
/// For T-shaped Pi-stacking, the distance threshold is set to 5.0 Å according to [getcontacts](https://getcontacts.github.io/interactions.html).
pub(super) fn find_pi_pi(ring1: &Plane, ring2: &Plane) -> Option<Interaction> {
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
                _ => {
                    if dist <= PI_T_DIST_THREHOLD {
                        Some(Interaction::PiTStacking)
                    } else {
                        None
                    }
                }
            },
            _ => None,
        }
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::contacts::residues::ResidueExt;
    use crate::load_model;
    use std::io::BufReader;

    #[test]
    fn unresolved_histidine_cation_pi_respects_protonation_policy() {
        let input = b"ATOM      1  CG  HIS A   1       0.000   0.000   3.000  1.00 20.00           C  \nEND\n";
        let pdb = ReadOptions::default()
            .set_format(Format::Pdb)
            .set_only_atomic_coords(true)
            .read_raw(BufReader::new(input.as_slice()))
            .unwrap()
            .0;
        let histidine = pdb.atoms_with_hierarchy().next().unwrap();
        let ring = Plane {
            center: na::Vector3::zeros(),
            normal: na::Vector3::new(0.0, 0.0, 1.0),
        };
        assert_eq!(
            find_cation_pi(&ring, &histidine, ProtonationMode::AllCharged, 7.4),
            Some(Interaction::PotentialCationPi)
        );
        assert_eq!(
            find_cation_pi(&ring, &histidine, ProtonationMode::ExplicitOnly, 7.4),
            None
        );
    }

    #[test]
    fn test_good_cation_pi() {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/6bft.pdb");
        let pdb = load_model(&path).unwrap().value;

        let y102_a = pdb
            .model(0)
            .unwrap()
            .chains()
            .find(|c| c.id() == "A")
            .unwrap()
            .residues()
            .find(|r| r.serial_number() == 102)
            .unwrap();
        let r82_g = pdb
            .atoms_with_hierarchy()
            .find(|x| {
                x.model().serial_number() == 0
                    && x.chain().id() == "G"
                    && x.conformer().name() == "ARG"
                    && x.residue().serial_number() == 82
                    && x.atom().name() == "NE"
            })
            .unwrap();
        let ring = y102_a.center_and_normal(Some(y102_a.ring_atoms())).unwrap();
        assert_eq!(
            find_cation_pi(&ring, &r82_g, ProtonationMode::AllCharged, 7.4),
            Some(Interaction::CationPi)
        );
    }

    #[test]
    fn test_bad_cation_pi_angle() {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/6bft.pdb");
        let pdb = load_model(&path).unwrap().value;

        let w108_a = pdb
            .model(0)
            .unwrap()
            .chains()
            .find(|c| c.id() == "A")
            .unwrap()
            .residues()
            .find(|r| r.serial_number() == 108)
            .unwrap();
        let k84_g = pdb
            .atoms_with_hierarchy()
            .find(|x| {
                x.model().serial_number() == 0
                    && x.chain().id() == "G"
                    && x.conformer().name() == "LYS"
                    && x.residue().serial_number() == 84
                    && x.atom().name() == "NZ"
            })
            .unwrap();
        let ring = w108_a.center_and_normal(Some(w108_a.ring_atoms())).unwrap();
        assert_eq!(
            find_cation_pi(&ring, &k84_g, ProtonationMode::AllCharged, 7.4),
            None
        );
    }
}
