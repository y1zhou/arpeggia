//! SC Calculator wrapping `SurfaceGenerator` with trimming and SC score calculation.
//! Ported from <https://github.com/cytokineking/sc-rs>

use core::f64;
use std::collections::{HashMap, HashSet};

use super::surface_generator::{SurfaceCalculatorError, SurfaceGenerator};
use super::types::*;
use super::vector3::{DotPoint, Vec3};
use super::{DOT_DENSITY, GAUSSIAN_WEIGHT, PERIPHERAL_BAND, SEPARATION_CUTOFF};
use pdbtbx::*;
use rayon::prelude::*;
use rstar::{PointDistance, RTree};

/// Clamp bounds for dot product to avoid numerical issues at extremes
const DOT_PRODUCT_CLAMP_MIN: f64 = -0.999;
const DOT_PRODUCT_CLAMP_MAX: f64 = 0.999;

#[derive(Default)]
pub struct ScCalculator {
    pub base: SurfaceGenerator,
}

impl ScCalculator {
    pub fn add_atoms(
        &mut self,
        pdb: &PDB,
        group1_chains: &HashSet<String>,
        group2_chains: &HashSet<String>,
    ) -> Result<(), SurfaceCalculatorError> {
        let tree = pdb.create_hierarchy_rtree();
        let max_radius_sq = SEPARATION_CUTOFF * SEPARATION_CUTOFF;

        let atoms = pdb
            .atoms_with_hierarchy()
            .filter_map(|hier| {
                let chain_id = hier.chain().id().to_string();

                // Determine which molecule (0 or 1) based on chain group
                let molecule = if group1_chains.contains(&chain_id) {
                    0
                } else if group2_chains.contains(&chain_id) {
                    1
                } else {
                    return None; // Skip atoms not in either group
                };

                let atom = hier.atom();
                let residue = hier.residue();
                let pos = atom.pos();

                let atom_radius = match SurfaceGenerator::get_atom_radius(
                    residue.name().unwrap_or(""),
                    atom.name(),
                ) {
                    Ok(r) => r,
                    // Fallback to pdbtbx element VdW radius
                    Err(_e) => {
                        match atom
                            .element()
                            .and_then(|element| element.atomic_radius().van_der_waals)
                        {
                            Some(r) => r,
                            None => {
                                return Some(Err(SurfaceCalculatorError::MissingRadius(format!(
                                    "atom {} ({}) has no usable van der Waals radius",
                                    atom.serial_number(),
                                    atom.name()
                                ))));
                            }
                        }
                    }
                };

                // Locate neighboring atoms within cutoff to assign attention
                // The items are (molecule index, Atom, distance squared)
                let neighbors: Vec<(usize, &Atom, f64)> = tree
                    .locate_within_distance(pos, max_radius_sq)
                    .map(|h| {
                        let h_chain = h.chain().id();
                        let h_mol = if group1_chains.contains(h_chain) {
                            0
                        } else if group2_chains.contains(h_chain) {
                            1
                        } else {
                            2 // Outside both groups, ignored later
                        };
                        (h_mol, h.atom(), atom.distance_2(&h.atom().pos()))
                    })
                    .collect();

                let attention = match neighbors
                    .iter()
                    .map(|(m, _, d)| {
                        if (m == &2) | (&molecule == m) {
                            f64::INFINITY
                        } else {
                            *d
                        }
                    })
                    .fold(f64::INFINITY, f64::min)
                {
                    d if d < max_radius_sq => Attention::Buried,
                    _ => Attention::Far,
                };

                Some(Ok(ScAtom {
                    atomi: atom.serial_number(),
                    molecule,
                    radius: atom_radius,
                    density: DOT_DENSITY,
                    attention,
                    accessible: false,
                    atomn: atom.name().to_string(),
                    resn: residue.name().unwrap_or("UNK").to_string(),
                    coor: Vec3::new(pos.0, pos.1, pos.2),
                    neighbors_atomi_dist2: neighbors
                        .into_iter()
                        .map(|(_, atom, d)| (atom.serial_number(), d))
                        .collect::<HashMap<usize, f64>>(),
                    neighbor_indices: Vec::new(),
                    buried_by_indices: Vec::new(),
                }))
            })
            .collect::<Result<Vec<_>, _>>()?;

        // Add atoms to calculator
        self.base.run.atoms.clear();
        self.base.run.atoms = atoms;
        Ok(())
    }

    pub fn calc(&mut self) -> Result<Results, SurfaceCalculatorError> {
        if self.base.run.atoms.is_empty() {
            return Err(SurfaceCalculatorError::NoAtoms);
        }
        for i in 0..2 {
            if !self.base.run.atoms.iter().any(|atom| atom.molecule == i) {
                return Err(SurfaceCalculatorError::NoAtomsForGroup(i + 1));
            }
        }

        self.base.generate_molecular_surfaces()?;

        if self.base.run.dots.iter().any(|x| x.is_empty()) {
            return Err(SurfaceCalculatorError::NoInterface);
        }

        for i in 0..2 {
            self.trim_peripheral_band(i);
        }
        if self.base.run.trimmed_dots.iter().any(Vec::is_empty) {
            return Err(SurfaceCalculatorError::NoInterface);
        }

        self.calc_neighbor_score(0, 1);
        self.calc_neighbor_score(1, 0);

        self.base.run.results.sc = f64::midpoint(
            self.base.run.results.surfaces[0].s_median,
            self.base.run.results.surfaces[1].s_median,
        );

        Ok(self.base.run.results.clone())
    }

    /// Trims peripheral dots that are within distance of non-buried dots.
    /// Uses RTree spatial indexing for efficient range queries.
    fn trim_peripheral_band(&mut self, i: usize) {
        let sdots = &self.base.run.dots[i];
        let r2 = PERIPHERAL_BAND * PERIPHERAL_BAND;

        // Build RTree of non-buried dots for efficient range queries
        let non_buried_points: Vec<DotPoint> = sdots
            .iter()
            .enumerate()
            .filter(|(_, dot)| !dot.buried)
            .map(|(idx, dot)| DotPoint::new(idx, dot.coor))
            .collect();

        let non_buried_tree: RTree<DotPoint> = RTree::bulk_load(non_buried_points);

        // Filter buried dots that are NOT within r of any non-buried dot
        let indices: Vec<usize> = (0..sdots.len())
            .into_par_iter()
            .filter(|&idx| {
                let dot = &sdots[idx];
                if !dot.buried {
                    return false;
                }
                let query_point = [dot.coor.x, dot.coor.y, dot.coor.z];
                // Check if there's any non-buried dot within distance r
                // Using locate_within_distance for efficient O(log n) lookup
                non_buried_tree
                    .locate_within_distance(query_point, r2)
                    .next()
                    .is_none()
            })
            .collect();

        self.base.run.trimmed_dots[i].clear();
        self.base.run.trimmed_dots[i] = indices;
    }

    /// Calculate normal-complementarity scores between two trimmed surfaces.
    /// Uses RTree spatial indexing for efficient nearest neighbor search.
    fn calc_neighbor_score(&mut self, my: usize, their: usize) {
        let my_dots = &self.base.run.trimmed_dots[my];
        let their_dots = &self.base.run.trimmed_dots[their];
        if my_dots.is_empty() || their_dots.is_empty() {
            return;
        }

        // Build RTree for buried dots on the opposite surface
        let their_buried_points: Vec<DotPoint> = their_dots
            .iter()
            .filter_map(|&idx| {
                let dot = &self.base.run.dots[their][idx];
                if dot.buried {
                    Some(DotPoint::new(idx, dot.coor))
                } else {
                    None
                }
            })
            .collect();

        if their_buried_points.is_empty() {
            return;
        }

        let their_tree: RTree<DotPoint> = RTree::bulk_load(their_buried_points);

        let run = &self.base.run;
        let mut scores: Vec<f64> = my_dots
            .par_iter()
            .filter_map(|&index| {
                let dot = &run.dots[my][index];
                let query = [dot.coor.x, dot.coor.y, dot.coor.z];
                their_tree.nearest_neighbor(&query).map(|nearest| {
                    let neighbor = &run.dots[their][nearest.index];
                    let distance_squared = nearest.distance_2(&query);
                    let score = dot.outnml.dot(neighbor.outnml)
                        * (-distance_squared * GAUSSIAN_WEIGHT).exp();
                    -score.clamp(DOT_PRODUCT_CLAMP_MIN, DOT_PRODUCT_CLAMP_MAX)
                })
            })
            .collect();
        if scores.is_empty() {
            return;
        }
        let median_index = scores.len() / 2;
        let (_, median, _) = scores.select_nth_unstable_by(median_index, |a, b| {
            a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
        });
        self.base.run.results.surfaces[my].s_median = *median;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::load_model;

    #[test]
    fn records_atoms_for_both_surface_groups() {
        let path = format!("{}/test-data/6bft.pdb", env!("CARGO_MANIFEST_DIR"));
        let pdb = load_model(&path).unwrap().value;
        let group1 = HashSet::from(["H".to_string()]);
        let group2 = HashSet::from(["L".to_string()]);
        let mut calculator = ScCalculator::default();

        calculator.add_atoms(&pdb, &group1, &group2).unwrap();

        assert!(
            calculator
                .base
                .run
                .atoms
                .iter()
                .any(|atom| atom.molecule == 0)
        );
        assert!(
            calculator
                .base
                .run
                .atoms
                .iter()
                .any(|atom| atom.molecule == 1)
        );
    }
}
