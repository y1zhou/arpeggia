//! SC Calculator wrapping `SurfaceGenerator` with trimming and SC score calculation.
//! Ported from <https://github.com/cytokineking/sc-rs>

use core::f64;
use std::collections::HashSet;

use super::surface_generator::{SurfaceCalculatorError, SurfaceGenerator};
use super::types::*;
use super::vector3::{DotPoint, Vec3};
use super::{GAUSSIAN_WEIGHT, PERIPHERAL_BAND, PROBE_RADIUS, SEPARATION_CUTOFF};
use pdbtbx::*;
use rayon::prelude::*;
use rstar::{PointDistance, RTree, primitives::GeomWithData};

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
        let max_radius_sq = SEPARATION_CUTOFF * SEPARATION_CUTOFF;

        let mut atoms = pdb
            .atoms_with_hierarchy()
            .filter_map(|hier| {
                let chain_id = hier.chain().id();

                // Determine which molecule (0 or 1) based on chain group
                let molecule = if group1_chains.contains(chain_id) {
                    0
                } else if group2_chains.contains(chain_id) {
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

                Some(Ok(ScAtom {
                    atomi: atom.serial_number(),
                    molecule,
                    radius: atom_radius,
                    attention: Attention::Far,
                    accessible: false,
                    atomn: atom.name().to_string(),
                    resn: residue.name().unwrap_or("UNK").to_string(),
                    coor: Vec3::new(pos.0, pos.1, pos.2),
                    neighbor_indices: Vec::new(),
                    neighbor_distances: Vec::new(),
                }))
            })
            .collect::<Result<Vec<_>, _>>()?;

        let tree = RTree::bulk_load(
            atoms
                .iter()
                .enumerate()
                .map(|(index, atom)| {
                    GeomWithData::new([atom.coor.x, atom.coor.y, atom.coor.z], index)
                })
                .collect(),
        );
        let neighbors = atoms
            .par_iter()
            .enumerate()
            .map(|(index, atom)| {
                let mut distances = Vec::new();
                let mut attention = Attention::Far;
                for neighbor in tree
                    .locate_within_distance([atom.coor.x, atom.coor.y, atom.coor.z], max_radius_sq)
                {
                    let neighbor_index = neighbor.data;
                    if neighbor_index == index {
                        continue;
                    }
                    let other = &atoms[neighbor_index];
                    let distance_squared = atom.coor.distance_squared(other.coor);
                    if atom.molecule != other.molecule {
                        attention = Attention::Buried;
                        continue;
                    }
                    if distance_squared <= 0.0001 {
                        return Err(SurfaceCalculatorError::Coincident(format!(
                            "{}:{}:{} == {}:{}:{}",
                            atom.atomi, atom.resn, atom.atomn, other.atomi, other.resn, other.atomn
                        )));
                    }
                    let bridge = atom.radius + other.radius + 2.0 * PROBE_RADIUS;
                    if distance_squared < bridge * bridge {
                        distances.push((neighbor_index, distance_squared));
                    }
                }
                let mut neighbor_indices = distances.clone();
                neighbor_indices.sort_unstable_by(|left, right| left.1.total_cmp(&right.1));
                distances.sort_unstable_by_key(|&(neighbor, _)| neighbor);
                Ok((
                    neighbor_indices
                        .into_iter()
                        .map(|(neighbor, _)| neighbor)
                        .collect::<Vec<_>>(),
                    distances,
                    attention,
                ))
            })
            .collect::<Result<Vec<_>, SurfaceCalculatorError>>()?;
        for (atom, (neighbor_indices, neighbor_distances, attention)) in
            atoms.iter_mut().zip(neighbors)
        {
            atom.accessible = neighbor_indices.is_empty();
            atom.neighbor_indices = neighbor_indices;
            atom.neighbor_distances = neighbor_distances;
            atom.attention = attention;
        }

        // Add atoms to calculator
        self.base.run.atoms = atoms;
        Ok(())
    }

    pub fn calc(&mut self) -> Result<Results, SurfaceCalculatorError> {
        if self.base.run.atoms.is_empty() {
            return Err(SurfaceCalculatorError::NoAtoms);
        }
        for i in 0..2 {
            if !self
                .base
                .run
                .atoms
                .iter()
                .any(|atom| usize::from(atom.molecule) == i)
            {
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
