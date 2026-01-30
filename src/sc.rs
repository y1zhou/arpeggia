//! Shape Complementarity (SC) calculations.
//!
//! This module provides functions for calculating shape complementarity (SC) between
//! two protein surfaces at their interface. The algorithm follows Lawrence & Colman (1993):
//! "Shape Complementarity at Protein/Protein Interfaces" (J Mol Biol 234:946-950).
//!
//! The SC score ranges from 0 to 1, where:
//! - 0 indicates no complementarity (surfaces don't match)
//! - 1 indicates perfect complementarity (surfaces fit perfectly together)
//!
//! Typical values for protein-protein interfaces range from 0.5 to 0.7.

use crate::sasa::{filter_pdb_by_model, prepare_pdb_for_sasa};
use crate::utils::parse_groups;
use nalgebra::Vector3;
use pdbtbx::*;
use rayon::prelude::*;
use std::collections::HashSet;

/// Default probe radius in Ångströms for molecular surface calculation
const DEFAULT_PROBE_RADIUS: f64 = 1.7;

/// Default dot density (dots per Å²) for surface generation
const DEFAULT_DOT_DENSITY: f64 = 15.0;

/// Default peripheral band width in Ångströms for trimming
const DEFAULT_BAND_WIDTH: f64 = 1.5;

/// Default separation distance in Ångströms for attention number assignment
const DEFAULT_SEPARATION: f64 = 8.0;

/// Weight factor for distance-weighted normal vector dot product
const DEFAULT_DISTANCE_WEIGHT: f64 = 0.5;

/// Result of a shape complementarity calculation.
#[derive(Debug, Clone)]
pub struct ScResult {
    /// The overall shape complementarity score (median, 0-1 scale)
    pub sc: f64,
    /// Mean distance between facing surfaces in Ångströms
    pub distance: f64,
    /// Total interface area (sum of both trimmed surfaces) in Ų
    pub area: f64,
    /// Per-surface results
    pub surfaces: [ScSurfaceResult; 2],
    /// Whether the calculation completed successfully
    pub valid: bool,
}

/// Result for a single molecular surface.
#[derive(Debug, Clone, Default)]
pub struct ScSurfaceResult {
    /// Number of atoms in this surface
    pub n_atoms: usize,
    /// Number of buried (interface) atoms
    pub n_buried_atoms: usize,
    /// Number of blocked (non-interface) atoms
    pub n_blocked_atoms: usize,
    /// Total number of surface dots
    pub n_dots: usize,
    /// Number of trimmed (interface) dots
    pub n_trimmed_dots: usize,
    /// Trimmed interface area in Ų
    pub trimmed_area: f64,
    /// Mean distance to other surface
    pub d_mean: f64,
    /// Median distance to other surface
    pub d_median: f64,
    /// Mean shape complementarity
    pub s_mean: f64,
    /// Median shape complementarity
    pub s_median: f64,
}

/// Configuration settings for SC calculation.
#[derive(Debug, Clone)]
pub struct ScSettings {
    /// Probe radius for surface generation (Å)
    pub probe_radius: f64,
    /// Dot density for surface (dots/Ų)
    pub density: f64,
    /// Peripheral band width for trimming (Å)
    pub band: f64,
    /// Separation distance for attention number assignment (Å)
    pub separation: f64,
    /// Weight for distance-weighted dot products
    pub weight: f64,
}

impl Default for ScSettings {
    fn default() -> Self {
        Self {
            probe_radius: DEFAULT_PROBE_RADIUS,
            density: DEFAULT_DOT_DENSITY,
            band: DEFAULT_BAND_WIDTH,
            separation: DEFAULT_SEPARATION,
            weight: DEFAULT_DISTANCE_WEIGHT,
        }
    }
}

/// A surface dot with position and normal vector.
#[derive(Debug, Clone)]
struct SurfaceDot {
    /// Position of the dot
    position: Vector3<f64>,
    /// Outward normal vector at this point
    normal: Vector3<f64>,
    /// Surface area represented by this dot
    area: f64,
    /// Which molecule this dot belongs to (0 or 1)
    molecule: usize,
    /// Whether this dot is buried (at interface)
    buried: bool,
}

/// An atom for SC calculation with radius and attention number.
#[derive(Debug, Clone)]
struct ScAtom {
    position: Vector3<f64>,
    radius: f64,
    molecule: usize,
    /// Whether this atom is at the interface (true) or blocked (false)
    is_interface: bool,
}

/// Generate surface dots for a collection of atoms using simplified Connolly surface.
fn generate_surface_dots(atoms: &[ScAtom], probe_radius: f64, density: f64) -> Vec<SurfaceDot> {
    atoms
        .par_iter()
        .enumerate()
        .flat_map(|(atom_idx, atom)| {
            if !atom.is_interface {
                return Vec::new();
            }

            let extended_radius = atom.radius + probe_radius;
            let surface_area = 4.0 * std::f64::consts::PI * extended_radius * extended_radius;
            let n_dots = (surface_area * density).ceil() as usize;
            let area_per_dot = surface_area / n_dots as f64;

            // Generate evenly distributed points on sphere using Fibonacci lattice
            let golden_ratio = (1.0 + 5.0_f64.sqrt()) / 2.0;
            let mut dots = Vec::with_capacity(n_dots);

            for i in 0..n_dots {
                let theta = 2.0 * std::f64::consts::PI * (i as f64) / golden_ratio;
                let phi = (1.0 - 2.0 * (i as f64 + 0.5) / n_dots as f64).acos();

                let normal = Vector3::new(
                    phi.sin() * theta.cos(),
                    phi.sin() * theta.sin(),
                    phi.cos(),
                );
                let position = atom.position + normal * extended_radius;

                // Check if this dot is occluded by any other atom
                let mut occluded = false;
                for (other_idx, other_atom) in atoms.iter().enumerate() {
                    if other_idx == atom_idx {
                        continue;
                    }
                    let other_extended = other_atom.radius + probe_radius;
                    let dist_sq = (position - other_atom.position).norm_squared();
                    if dist_sq < other_extended * other_extended {
                        occluded = true;
                        break;
                    }
                }

                if !occluded {
                    dots.push(SurfaceDot {
                        position,
                        normal,
                        area: area_per_dot,
                        molecule: atom.molecule,
                        buried: false, // Will be set later
                    });
                }
            }
            dots
        })
        .collect()
}

/// Classify dots as buried (interface) or exposed based on proximity to other molecule.
fn classify_dots_buried(dots: &mut [SurfaceDot], atoms: &[ScAtom], probe_radius: f64) {
    let mol0_positions: Vec<(Vector3<f64>, f64)> = atoms
        .iter()
        .filter(|a| a.molecule == 0)
        .map(|a| (a.position, a.radius))
        .collect();
    let mol1_positions: Vec<(Vector3<f64>, f64)> = atoms
        .iter()
        .filter(|a| a.molecule == 1)
        .map(|a| (a.position, a.radius))
        .collect();

    dots.par_iter_mut().for_each(|dot| {
        let other_atoms = if dot.molecule == 0 {
            &mol1_positions
        } else {
            &mol0_positions
        };

        // Check if any atom from the other molecule is close to this dot
        for (other_pos, other_radius) in other_atoms.iter() {
            let extended_radius = other_radius + probe_radius;
            let dist_sq = (dot.position - other_pos).norm_squared();
            // Use a looser criterion for "buried" - within 1.5x the extended radius
            if dist_sq < extended_radius * extended_radius * 2.25 {
                dot.buried = true;
                break;
            }
        }
    });
}

/// Trim peripheral band: keep only buried dots that have accessible dots nearby.
fn trim_peripheral_band(dots: &[SurfaceDot], band_width: f64) -> Vec<&SurfaceDot> {
    let band_sq = band_width * band_width;

    // Separate buried and accessible dots for this molecule
    let buried_dots: Vec<_> = dots.iter().filter(|d| d.buried).collect();
    let accessible_dots: Vec<_> = dots.iter().filter(|d| !d.buried).collect();

    if accessible_dots.is_empty() {
        // If no accessible dots, all buried dots are in the interface
        return buried_dots;
    }

    // Keep buried dots that have an accessible neighbor within band distance
    buried_dots
        .into_par_iter()
        .filter(|dot| {
            for acc_dot in &accessible_dots {
                let dist_sq = (dot.position - acc_dot.position).norm_squared();
                if dist_sq <= band_sq {
                    return true;
                }
            }
            false
        })
        .collect()
}

/// Calculate shape complementarity statistics between two surfaces.
fn calculate_sc_statistics(
    my_dots: &[&SurfaceDot],
    their_dots: &[&SurfaceDot],
    weight: f64,
) -> (f64, f64, f64, f64) {
    if my_dots.is_empty() || their_dots.is_empty() {
        return (0.0, 0.0, 0.0, 0.0);
    }

    let results: Vec<(f64, f64)> = my_dots
        .par_iter()
        .map(|my_dot| {
            // Find closest neighbor in their surface
            let mut min_dist_sq = f64::MAX;
            let mut closest_idx = 0;

            for (i, their_dot) in their_dots.iter().enumerate() {
                let dist_sq = (my_dot.position - their_dot.position).norm_squared();
                if dist_sq < min_dist_sq {
                    min_dist_sq = dist_sq;
                    closest_idx = i;
                }
            }

            let dist = min_dist_sq.sqrt();

            // Calculate dot product of normals, weighted by distance
            let their_dot = their_dots[closest_idx];
            let dot_product = my_dot.normal.dot(&their_dot.normal);

            // Weight by exponential of distance squared
            let weighted = dot_product * (-dist * dist * weight).exp();
            let clamped = weighted.clamp(-0.999, 0.999);

            (dist, clamped)
        })
        .collect();

    let distances: Vec<f64> = results.iter().map(|(d, _)| *d).collect();
    let weighted_products: Vec<f64> = results.iter().map(|(_, w)| *w).collect();

    // Calculate statistics
    let n = distances.len() as f64;
    let d_mean = distances.iter().sum::<f64>() / n;
    let s_mean = -weighted_products.iter().sum::<f64>() / n;

    // Sort for median calculation
    let mut sorted_distances = distances.clone();
    let mut sorted_products = weighted_products.clone();
    sorted_distances.sort_by(|a, b| a.partial_cmp(b).unwrap());
    sorted_products.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let d_median = if sorted_distances.len() % 2 == 0 {
        let mid = sorted_distances.len() / 2;
        (sorted_distances[mid - 1] + sorted_distances[mid]) / 2.0
    } else {
        sorted_distances[sorted_distances.len() / 2]
    };

    let s_median = if sorted_products.len() % 2 == 0 {
        let mid = sorted_products.len() / 2;
        -(sorted_products[mid - 1] + sorted_products[mid]) / 2.0
    } else {
        -sorted_products[sorted_products.len() / 2]
    };

    (d_mean, d_median, s_mean, s_median)
}

/// Calculate shape complementarity between two chain groups.
///
/// This function calculates the Lawrence & Colman shape complementarity score
/// between two sets of protein chains at their interface.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `groups` - Chain groups specification (e.g., "A,B/C,D")
/// * `model_num` - Model number to analyze (0 for first model)
/// * `settings` - Optional SC settings (uses defaults if None)
///
/// # Returns
///
/// An `ScResult` struct containing the SC score and related statistics.
///
/// # Example
///
/// ```no_run
/// use arpeggia::{load_model, get_sc};
///
/// let input_file = "path/to/structure.pdb".to_string();
/// let (pdb, _errors) = load_model(&input_file);
///
/// // Calculate SC between chains A,B and C,D
/// let result = get_sc(&pdb, "A,B/C,D", 0, None);
/// println!("Shape complementarity: {:.3}", result.sc);
/// ```
pub fn get_sc(pdb: &PDB, groups: &str, model_num: usize, settings: Option<ScSettings>) -> ScResult {
    let settings = settings.unwrap_or_default();

    // Prepare PDB: remove solvent, ions, hydrogens
    let pdb_prepared = prepare_pdb_for_sasa(pdb, true, true, "");
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    // Parse chain groups using the utility function
    let all_chains: HashSet<String> = pdb_filtered.chains().map(|c| c.id().to_string()).collect();
    let (group1_chains, group2_chains) = parse_groups(&all_chains, groups);

    // Use pdbtbx's RTree for efficient spatial queries
    let tree = pdb_filtered.create_hierarchy_rtree();
    let separation_sq = settings.separation * settings.separation;

    // Collect atoms with their molecule assignments and interface status
    let mut atoms: Vec<ScAtom> = Vec::new();
    let mut surface_results = [ScSurfaceResult::default(), ScSurfaceResult::default()];

    for hier in pdb_filtered.atoms_with_hierarchy() {
        let chain_id = hier.chain().id().to_string();
        let molecule = if group1_chains.contains(&chain_id) {
            0
        } else if group2_chains.contains(&chain_id) {
            1
        } else {
            continue; // Skip chains not in either group
        };

        // Skip hydrogen atoms
        if hier.atom().element() == Some(&Element::H) {
            continue;
        }

        let vdw_radius = hier
            .atom()
            .element()
            .and_then(|e| e.atomic_radius().van_der_waals)
            .unwrap_or(1.5);

        let pos = hier.atom().pos();
        let position = Vector3::new(pos.0, pos.1, pos.2);

        // Check if this atom is at the interface by looking for nearby atoms from the other molecule
        let other_group = if molecule == 0 {
            &group2_chains
        } else {
            &group1_chains
        };

        let is_interface = tree
            .locate_within_distance(pos, separation_sq)
            .any(|neighbor| other_group.contains(neighbor.chain().id()));

        surface_results[molecule].n_atoms += 1;
        if is_interface {
            surface_results[molecule].n_buried_atoms += 1;
        } else {
            surface_results[molecule].n_blocked_atoms += 1;
        }

        atoms.push(ScAtom {
            position,
            radius: vdw_radius,
            molecule,
            is_interface,
        });
    }

    if atoms.is_empty() || surface_results[0].n_atoms == 0 || surface_results[1].n_atoms == 0 {
        return ScResult {
            sc: 0.0,
            distance: 0.0,
            area: 0.0,
            surfaces: surface_results,
            valid: false,
        };
    }

    // Generate molecular surface dots
    let mut dots = generate_surface_dots(&atoms, settings.probe_radius, settings.density);

    if dots.is_empty() {
        return ScResult {
            sc: 0.0,
            distance: 0.0,
            area: 0.0,
            surfaces: surface_results,
            valid: false,
        };
    }

    // Classify dots as buried (interface) or exposed
    classify_dots_buried(&mut dots, &atoms, settings.probe_radius);

    // Separate dots by molecule
    let mol0_dots: Vec<SurfaceDot> = dots.iter().filter(|d| d.molecule == 0).cloned().collect();
    let mol1_dots: Vec<SurfaceDot> = dots.iter().filter(|d| d.molecule == 1).cloned().collect();

    surface_results[0].n_dots = mol0_dots.len();
    surface_results[1].n_dots = mol1_dots.len();

    // Trim peripheral band for each surface
    let trimmed_mol0 = trim_peripheral_band(&mol0_dots, settings.band);
    let trimmed_mol1 = trim_peripheral_band(&mol1_dots, settings.band);

    surface_results[0].n_trimmed_dots = trimmed_mol0.len();
    surface_results[1].n_trimmed_dots = trimmed_mol1.len();
    surface_results[0].trimmed_area = trimmed_mol0.iter().map(|d| d.area).sum();
    surface_results[1].trimmed_area = trimmed_mol1.iter().map(|d| d.area).sum();

    if trimmed_mol0.is_empty() || trimmed_mol1.is_empty() {
        return ScResult {
            sc: 0.0,
            distance: 0.0,
            area: 0.0,
            surfaces: surface_results,
            valid: false,
        };
    }

    // Calculate statistics for each direction
    let (d_mean_0, d_median_0, s_mean_0, s_median_0) =
        calculate_sc_statistics(&trimmed_mol0, &trimmed_mol1, settings.weight);
    let (d_mean_1, d_median_1, s_mean_1, s_median_1) =
        calculate_sc_statistics(&trimmed_mol1, &trimmed_mol0, settings.weight);

    surface_results[0].d_mean = d_mean_0;
    surface_results[0].d_median = d_median_0;
    surface_results[0].s_mean = s_mean_0;
    surface_results[0].s_median = s_median_0;

    surface_results[1].d_mean = d_mean_1;
    surface_results[1].d_median = d_median_1;
    surface_results[1].s_mean = s_mean_1;
    surface_results[1].s_median = s_median_1;

    // Average the results from both directions
    let sc = (s_median_0 + s_median_1) / 2.0;
    let distance = (d_median_0 + d_median_1) / 2.0;
    let area = surface_results[0].trimmed_area + surface_results[1].trimmed_area;

    ScResult {
        sc,
        distance,
        area,
        surfaces: surface_results,
        valid: true,
    }
}

/// Get the shape complementarity score as a simple f64 value.
///
/// This is a convenience wrapper around `get_sc` that returns just the SC score.
///
/// # Arguments
///
/// * `pdb` - Reference to a PDB structure
/// * `groups` - Chain groups specification (e.g., "A,B/C,D")
/// * `model_num` - Model number to analyze (0 for first model)
///
/// # Returns
///
/// The SC score (0-1), or -1.0 if the calculation failed.
pub fn get_sc_score(pdb: &PDB, groups: &str, model_num: usize) -> f64 {
    let result = get_sc(pdb, groups, model_num, None);
    if result.valid {
        result.sc
    } else {
        -1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::load_model;

    fn load_multi_chain() -> PDB {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/6bft.pdb");
        let (pdb, _) = load_model(&path);
        pdb
    }

    #[test]
    fn test_sc_settings_default() {
        let settings = ScSettings::default();
        assert!((settings.probe_radius - 1.7).abs() < 0.001);
        assert!((settings.density - 15.0).abs() < 0.001);
        assert!((settings.band - 1.5).abs() < 0.001);
        assert!((settings.separation - 8.0).abs() < 0.001);
    }

    #[test]
    fn test_sc_multi_chain() {
        let pdb = load_multi_chain();

        // 6bft has chains: A, B, C, G, H, L (antibody-antigen complex)
        let result = get_sc(&pdb, "H,L/A,B,C", 0, None);

        // Should complete successfully
        assert!(result.valid, "SC calculation should be valid");

        // SC score should be between 0 and 1
        assert!(
            result.sc >= 0.0 && result.sc <= 1.0,
            "SC score should be between 0 and 1, got {}",
            result.sc
        );

        // Should have some interface area
        assert!(result.area > 0.0, "Should have positive interface area");

        // Both surfaces should have some trimmed dots
        assert!(
            result.surfaces[0].n_trimmed_dots > 0,
            "Surface 0 should have trimmed dots"
        );
        assert!(
            result.surfaces[1].n_trimmed_dots > 0,
            "Surface 1 should have trimmed dots"
        );
    }

    #[test]
    fn test_sc_score_wrapper() {
        let pdb = load_multi_chain();
        let score = get_sc_score(&pdb, "H,L/A,B,C", 0);

        // Should return valid score
        assert!(score >= 0.0, "SC score should be >= 0");
        assert!(score <= 1.0, "SC score should be <= 1");
    }

    #[test]
    fn test_surface_dot_generation() {
        let atoms = vec![ScAtom {
            position: Vector3::new(0.0, 0.0, 0.0),
            radius: 1.5,
            molecule: 0,
            is_interface: true,
        }];

        let dots = generate_surface_dots(&atoms, 1.7, 15.0);

        // Should generate some dots
        assert!(!dots.is_empty(), "Should generate surface dots");

        // All dots should have the same molecule assignment
        assert!(dots.iter().all(|d| d.molecule == 0));

        // Dots should be roughly at the expected distance from center
        let expected_dist = 1.5 + 1.7; // atom radius + probe radius
        for dot in &dots {
            let dist = dot.position.norm();
            assert!(
                (dist - expected_dist).abs() < 0.01,
                "Dot distance {} should be close to {}",
                dist,
                expected_dist
            );
        }
    }
}
