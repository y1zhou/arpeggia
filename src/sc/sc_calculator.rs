//! SC Calculator wrapping SurfaceGenerator with trimming and SC score calculation.
//! Ported from https://github.com/cytokineking/sc-rs

use super::surface_generator::{SurfaceCalculatorError, SurfaceGenerator};
use super::types::*;
use super::vector3::ScValue;
use rayon::prelude::*;

/// Large initial distance squared for neighbor search
const MAX_DISTANCE_SQUARED: f64 = 9.0e20;
/// Clamp bounds for dot product to avoid numerical issues at extremes
const DOT_PRODUCT_CLAMP_MIN: f64 = -0.999;
const DOT_PRODUCT_CLAMP_MAX: f64 = 0.999;

pub struct ScCalculator {
    pub base: SurfaceGenerator,
}

impl Default for ScCalculator {
    fn default() -> Self {
        Self::new()
    }
}

impl ScCalculator {
    pub fn new() -> Self {
        Self {
            base: SurfaceGenerator::new(),
        }
    }

    /// Set the number of threads for parallel calculations.
    /// 0 means automatic (use all available cores).
    pub fn set_threads(&mut self, threads: usize) {
        self.base.settings.threads = threads;
        self.base.settings.enable_parallel = threads != 1;
    }

    pub fn calc(&mut self) -> Result<Results, SurfaceCalculatorError> {
        self.base.init()?;
        self.base.run.results.valid = 0;

        if self.base.run.atoms.is_empty() {
            return Err(SurfaceCalculatorError::NoAtoms);
        }
        if self.base.run.results.surfaces[0].n_atoms == 0 {
            return Err(SurfaceCalculatorError::Io(std::io::Error::other(
                "No atoms for molecule 1",
            )));
        }
        if self.base.run.results.surfaces[1].n_atoms == 0 {
            return Err(SurfaceCalculatorError::Io(std::io::Error::other(
                "No atoms for molecule 2",
            )));
        }

        self.base.assign_attention_numbers();
        self.base.generate_molecular_surfaces()?;

        if self.base.run.dots[0].is_empty() || self.base.run.dots[1].is_empty() {
            return Err(SurfaceCalculatorError::Io(std::io::Error::other(
                "No molecular dots generated",
            )));
        }

        for i in 0..2 {
            let area = self.trim_peripheral_band(i)?;
            self.base.run.results.surfaces[i].trimmed_area = area;
            self.base.run.results.surfaces[i].n_trimmed_dots = self.base.run.trimmed_dots[i].len();
            self.base.run.results.surfaces[i].n_all_dots = self.base.run.dots[i].len();
        }

        self.calc_neighbor_distance(0, 1);
        self.calc_neighbor_distance(1, 0);

        // Combine results from both surfaces
        self.base.run.results.combined.d_mean = (self.base.run.results.surfaces[0].d_mean
            + self.base.run.results.surfaces[1].d_mean)
            / 2.0;
        self.base.run.results.combined.d_median = (self.base.run.results.surfaces[0].d_median
            + self.base.run.results.surfaces[1].d_median)
            / 2.0;
        self.base.run.results.combined.s_mean = (self.base.run.results.surfaces[0].s_mean
            + self.base.run.results.surfaces[1].s_mean)
            / 2.0;
        self.base.run.results.combined.s_median = (self.base.run.results.surfaces[0].s_median
            + self.base.run.results.surfaces[1].s_median)
            / 2.0;
        self.base.run.results.combined.n_atoms = self.base.run.results.surfaces[0].n_atoms
            + self.base.run.results.surfaces[1].n_atoms;
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

    fn trim_peripheral_band(&mut self, i: usize) -> Result<ScValue, SurfaceCalculatorError> {
        let (indices, area) = if self.base.settings.enable_parallel {
            let sdots = &self.base.run.dots[i];
            self.base.run_parallel(|| {
                let indices: Vec<usize> = (0..sdots.len())
                    .into_par_iter()
                    .filter(|&idx| sdots[idx].buried && self.trim_peripheral_band_check_dot(idx, sdots))
                    .collect();
                let area: f64 = indices.par_iter().map(|&idx| sdots[idx].area).sum();
                (indices, area)
            })
        } else {
            let sdots = &self.base.run.dots[i];
            let indices: Vec<usize> = sdots
                .iter()
                .enumerate()
                .filter(|(idx, dot)| dot.buried && self.trim_peripheral_band_check_dot(*idx, sdots))
                .map(|(idx, _)| idx)
                .collect();
            let area: f64 = indices.iter().map(|&idx| sdots[idx].area).sum();
            (indices, area)
        };

        self.base.run.trimmed_dots[i].clear();
        self.base.run.trimmed_dots[i] = indices;
        Ok(area)
    }

    fn trim_peripheral_band_check_dot(&self, dot_index: usize, sdots: &[Dot]) -> bool {
        let r2 = self.base.settings.peripheral_band * self.base.settings.peripheral_band;
        let dot = &sdots[dot_index];
        !sdots.iter().enumerate().any(|(i, dot2)| {
            i != dot_index && !dot2.buried && dot.coor.distance_squared(dot2.coor) <= r2
        })
    }

    fn calc_neighbor_distance(&mut self, my: usize, their: usize) {
        let my_dots = &self.base.run.trimmed_dots[my];
        let their_dots = &self.base.run.trimmed_dots[their];
        if my_dots.is_empty() || their_dots.is_empty() {
            return;
        }

        let (distances, scores, distmin_sum, score_sum) = if self.base.settings.enable_parallel {
            let gaussian_w = self.base.settings.gaussian_w;
            let run_ref = &self.base.run;
            self.base.run_parallel(|| {
                let pairs: Vec<(f64, f64)> = my_dots
                    .par_iter()
                    .filter_map(|&pd| {
                        let dot1 = &run_ref.dots[my][pd];
                        let (distmin2, neighbor) = their_dots.iter().fold(
                            (MAX_DISTANCE_SQUARED, None),
                            |(min_d2, min_n), &pd2| {
                                let dot2 = &run_ref.dots[their][pd2];
                                if !dot2.buried {
                                    return (min_d2, min_n);
                                }
                                let d2 = dot2.coor.distance_squared(dot1.coor);
                                if d2 <= min_d2 {
                                    (d2, Some(dot2))
                                } else {
                                    (min_d2, min_n)
                                }
                            },
                        );
                        neighbor.map(|n| {
                            let distmin = distmin2.sqrt();
                            let mut r = dot1.outnml.dot(n.outnml);
                            r *= (-(distmin * distmin) * gaussian_w).exp();
                            r = r.clamp(DOT_PRODUCT_CLAMP_MIN, DOT_PRODUCT_CLAMP_MAX);
                            (distmin, -r)
                        })
                    })
                    .collect();
                let (distances, scores): (Vec<f64>, Vec<f64>) = pairs.into_iter().unzip();
                let distmin_sum: f64 = distances.par_iter().sum();
                let score_sum: f64 = scores.par_iter().map(|v| -v).sum();
                (distances, scores, distmin_sum, score_sum)
            })
        } else {
            let gaussian_w = self.base.settings.gaussian_w;
            let pairs: Vec<(f64, f64)> = my_dots
                .iter()
                .filter_map(|&pd| {
                    let dot1 = &self.base.run.dots[my][pd];
                    let (distmin2, neighbor) = their_dots.iter().fold(
                        (MAX_DISTANCE_SQUARED, None),
                        |(min_d2, min_n), &pd2| {
                            let dot2 = &self.base.run.dots[their][pd2];
                            if !dot2.buried {
                                return (min_d2, min_n);
                            }
                            let d2 = dot2.coor.distance_squared(dot1.coor);
                            if d2 <= min_d2 {
                                (d2, Some(dot2))
                            } else {
                                (min_d2, min_n)
                            }
                        },
                    );
                    neighbor.map(|n| {
                        let distmin = distmin2.sqrt();
                        let mut r = dot1.outnml.dot(n.outnml);
                        r *= (-(distmin * distmin) * gaussian_w).exp();
                        r = r.clamp(DOT_PRODUCT_CLAMP_MIN, DOT_PRODUCT_CLAMP_MAX);
                        (distmin, -r)
                    })
                })
                .collect();
            let (distances, scores): (Vec<f64>, Vec<f64>) = pairs.into_iter().unzip();
            let distmin_sum: f64 = distances.iter().sum();
            let score_sum: f64 = scores.iter().map(|v| -v).sum();
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

    pub fn add_atom(&mut self, molecule: i32, atom: Atom) -> Result<(), SurfaceCalculatorError> {
        self.base.add_atom(molecule, atom)
    }
}
