//! SC Calculator wrapping `SurfaceGenerator` with trimming and SC score calculation.
//! Ported from <https://github.com/cytokineking/sc-rs>

use core::f64;
use std::collections::{HashMap, HashSet};

use super::surface_generator::{SurfaceCalculatorError, SurfaceGenerator};
use super::types::*;
use super::vector3::{DotPoint, Vec3};
use pdbtbx::*;
use rayon::prelude::*;
use rstar::{PointDistance, RTree};

/// Clamp bounds for dot product to avoid numerical issues at extremes
const DOT_PRODUCT_CLAMP_MIN: f64 = -0.999;
const DOT_PRODUCT_CLAMP_MAX: f64 = 0.999;

pub struct ScCalculator {
    pub base: SurfaceGenerator,
}

impl ScCalculator {
    pub fn new() -> Self {
        Self {
            base: SurfaceGenerator::new(),
        }
    }

    pub fn add_atoms(
        &mut self,
        pdb: &PDB,
        group1_chains: &HashSet<String>,
        group2_chains: &HashSet<String>,
    ) -> Result<(), SurfaceCalculatorError> {
        self.base.init();

        let tree = pdb.create_hierarchy_rtree();
        let max_radius_sq =
            self.base.settings.separation_cutoff * self.base.settings.separation_cutoff;

        let atoms: Vec<ScAtom> = pdb
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

                let atom_radius = match self
                    .base
                    .get_atom_radius(residue.name().unwrap(), atom.name())
                {
                    Ok(r) => r,
                    // Fallback to pdbtbx element VdW radius
                    Err(_e) => {
                        match atom.element().unwrap().atomic_radius().van_der_waals {
                            Some(r) => r,
                            None => return None, // Skip if no radius info
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
                    .reduce(f64::min)
                    .unwrap()
                {
                    d if d < max_radius_sq => {
                        self.base.run.results.surfaces[molecule].n_buried_atoms += 1;
                        Attention::Buried
                    }
                    _ => {
                        self.base.run.results.surfaces[molecule].n_blocked_atoms += 1;
                        Attention::Far
                    }
                };

                Some(ScAtom {
                    atomi: atom.serial_number(),
                    molecule,
                    radius: atom_radius,
                    density: self.base.settings.dot_density,
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
                })
            })
            .collect();

        // Add atoms to calculator
        self.base.run.results.n_atoms = atoms.len();
        for mol in 0..1 {
            self.base.run.results.surfaces[mol].n_atoms =
                atoms.par_iter().filter(|atom| atom.molecule == mol).count();
        }
        self.base.run.atoms.clear();
        self.base.run.atoms = atoms;
        Ok(())
    }

    pub fn calc(&mut self) -> Result<Results, SurfaceCalculatorError> {
        self.base.init();
        self.base.run.results.valid = 0;

        if self.base.run.atoms.is_empty() {
            return Err(SurfaceCalculatorError::NoAtoms);
        }
        for i in 0..1 {
            if self.base.run.results.surfaces[i].n_atoms == 0 {
                return Err(SurfaceCalculatorError::Io(std::io::Error::other(format!(
                    "No atoms for chain group {}",
                    i + 1
                ))));
            }
        }

        self.base.generate_molecular_surfaces()?;

        if self.base.run.dots.iter().any(|x| x.is_empty()) {
            return Err(SurfaceCalculatorError::Io(std::io::Error::other(
                "No molecular dots generated",
            )));
        }

        for i in 0..2 {
            let area = self.trim_peripheral_band(i);
            self.base.run.results.surfaces[i].trimmed_area = area;
            self.base.run.results.surfaces[i].n_trimmed_dots = self.base.run.trimmed_dots[i].len();
            self.base.run.results.surfaces[i].n_all_dots = self.base.run.dots[i].len();
        }

        self.calc_neighbor_distance(0, 1);
        self.calc_neighbor_distance(1, 0);

        // Combine results from both surfaces
        self.base.run.results.combined.d_mean = f64::midpoint(
            self.base.run.results.surfaces[0].d_mean,
            self.base.run.results.surfaces[1].d_mean,
        );
        self.base.run.results.combined.d_median = f64::midpoint(
            self.base.run.results.surfaces[0].d_median,
            self.base.run.results.surfaces[1].d_median,
        );
        self.base.run.results.combined.s_mean = f64::midpoint(
            self.base.run.results.surfaces[0].s_mean,
            self.base.run.results.surfaces[1].s_mean,
        );
        self.base.run.results.combined.s_median = f64::midpoint(
            self.base.run.results.surfaces[0].s_median,
            self.base.run.results.surfaces[1].s_median,
        );
        self.base.run.results.combined.n_atoms =
            self.base.run.results.surfaces[0].n_atoms + self.base.run.results.surfaces[1].n_atoms;
        self.base.run.results.combined.n_buried_atoms = self.base.run.results.surfaces[0]
            .n_buried_atoms
            + self.base.run.results.surfaces[1].n_buried_atoms;
        self.base.run.results.combined.n_blocked_atoms = self.base.run.results.surfaces[0]
            .n_blocked_atoms
            + self.base.run.results.surfaces[1].n_blocked_atoms;
        self.base.run.results.combined.n_all_dots = self.base.run.results.surfaces[0].n_all_dots
            + self.base.run.results.surfaces[1].n_all_dots;
        self.base.run.results.combined.n_trimmed_dots = self.base.run.results.surfaces[0]
            .n_trimmed_dots
            + self.base.run.results.surfaces[1].n_trimmed_dots;
        self.base.run.results.combined.trimmed_area = self.base.run.results.surfaces[0]
            .trimmed_area
            + self.base.run.results.surfaces[1].trimmed_area;

        self.base.run.results.sc = self.base.run.results.combined.s_median;
        self.base.run.results.distance = self.base.run.results.combined.d_median;
        self.base.run.results.area = self.base.run.results.combined.trimmed_area;
        self.base.run.results.valid = 1;

        Ok(self.base.run.results.clone())
    }

    /// Trims peripheral dots that are within distance of non-buried dots.
    /// Uses RTree spatial indexing for efficient range queries.
    fn trim_peripheral_band(&mut self, i: usize) -> f64 {
        let sdots = &self.base.run.dots[i];
        let r2 = self.base.settings.peripheral_band * self.base.settings.peripheral_band;

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

        let area: f64 = indices.par_iter().map(|&idx| sdots[idx].area).sum();

        self.base.run.trimmed_dots[i].clear();
        self.base.run.trimmed_dots[i] = indices;
        area
    }

    /// Calculate distances and scores between trimmed dots from two surfaces.
    /// Uses RTree spatial indexing for efficient nearest neighbor search.
    fn calc_neighbor_distance(&mut self, my: usize, their: usize) {
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

        let (distances, scores, distmin_sum, score_sum) = {
            let gaussian_w = self.base.settings.gaussian_w;
            let run_ref = &self.base.run;
            let pairs: Vec<(f64, f64)> = my_dots
                .par_iter()
                .filter_map(|&pd| {
                    let dot1 = &run_ref.dots[my][pd];
                    let query_point = [dot1.coor.x, dot1.coor.y, dot1.coor.z];

                    // Use RTree nearest neighbor query for O(log n) lookup
                    their_tree
                        .nearest_neighbor(&query_point)
                        .map(|nearest_point| {
                            let neighbor = &run_ref.dots[their][nearest_point.index];
                            let distmin2 = nearest_point.distance_2(&query_point);
                            let distmin = distmin2.sqrt();
                            let mut r = dot1.outnml.dot(neighbor.outnml);
                            r *= (-distmin2 * gaussian_w).exp();
                            r = r.clamp(DOT_PRODUCT_CLAMP_MIN, DOT_PRODUCT_CLAMP_MAX);
                            (distmin, -r)
                        })
                })
                .collect();
            let (distances, scores): (Vec<f64>, Vec<f64>) = pairs.into_iter().unzip();
            let distmin_sum: f64 = distances.par_iter().sum();
            let score_sum: f64 = scores.par_iter().map(|v| -v).sum();
            (distances, scores, distmin_sum, score_sum)
        };

        let total_points = distances.len() as f64;
        if total_points == 0.0 {
            return;
        }

        let mut distances = distances;
        let mut scores = scores;
        let d_len = distances.len() as f64;
        let s_len = scores.len() as f64;

        // Calculate medians
        let d_median_val = {
            let median_idx = distances.len() / 2;
            let (_, m, _) = distances.select_nth_unstable_by(median_idx, |a, b| {
                a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
            });
            *m
        };
        let s_median_val = {
            let median_idx = scores.len() / 2;
            let (_, m, _) = scores.select_nth_unstable_by(median_idx, |a, b| {
                a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal)
            });
            *m
        };

        self.base.run.results.surfaces[my].d_mean = distmin_sum / d_len;
        self.base.run.results.surfaces[my].d_median = d_median_val;
        self.base.run.results.surfaces[my].s_mean = -(score_sum / s_len);
        self.base.run.results.surfaces[my].s_median = s_median_val;
    }
}
