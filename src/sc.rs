//! Shape Complementarity (SC) calculations.
//!
//! This module provides functions for calculating shape complementarity (SC) between
//! two protein surfaces at their interface. The algorithm follows Lawrence & Colman (1993):
//! "Shape Complementarity at Protein/Protein Interfaces" (J Mol Biol 234:946-950).

use crate::sasa::{filter_pdb_by_model, prepare_pdb_for_sasa};
use crate::utils::parse_groups;
use nalgebra::Vector3;
use pdbtbx::*;
use rayon::prelude::*;
use rstar::{RTree, RTreeObject, PointDistance, AABB};
use std::collections::HashSet;

/// Gaussian weight w = 0.5 Å^-2 (Lawrence & Colman 1993, Fig. 1)
const GAUSSIAN_W: f64 = 0.5;
/// Peripheral exclusion band d = 1.5 Å (Lawrence & Colman 1993)
const PERIPHERAL_BAND: f64 = 1.5;
/// Target dot density ~15 dots per Å^2 (Lawrence & Colman 1993)
const DOT_DENSITY: f64 = 15.0;
/// Probe radius (Connolly 1983)
const PROBE_RADIUS: f64 = 1.7;
/// Separation cutoff for interface classification
const SEPARATION_CUTOFF: f64 = 8.0;

/// Result of a shape complementarity calculation.
#[derive(Debug, Clone)]
pub struct ScResult {
    /// The overall shape complementarity score
    pub sc: f64,
    /// Median distance between facing surfaces
    pub distance: f64,
    /// Total interface area
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
    /// Trimmed interface area
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
    /// Gaussian weight for distance-weighted dot products (Å^-2)
    pub weight: f64,
}

impl Default for ScSettings {
    fn default() -> Self {
        Self {
            probe_radius: PROBE_RADIUS,
            density: DOT_DENSITY,
            band: PERIPHERAL_BAND,
            separation: SEPARATION_CUTOFF,
            weight: GAUSSIAN_W,
        }
    }
}

/// A surface dot with position and normal vector.
#[derive(Debug, Clone)]
struct Dot {
    /// Position of the dot (on VdW surface)
    coor: Vector3<f64>,
    /// Outward normal vector
    outnml: Vector3<f64>,
    /// Surface area represented by this dot
    area: f64,
    /// Whether this dot is buried (within opposite molecule's probe sphere)
    buried: bool,
}

/// Dot wrapper for RTree
#[derive(Debug, Clone)]
struct DotPoint {
    idx: usize,
    pos: [f64; 3],
}

impl RTreeObject for DotPoint {
    type Envelope = AABB<[f64; 3]>;
    fn envelope(&self) -> Self::Envelope {
        AABB::from_point(self.pos)
    }
}

impl PointDistance for DotPoint {
    fn distance_2(&self, point: &[f64; 3]) -> f64 {
        let dx = self.pos[0] - point[0];
        let dy = self.pos[1] - point[1];
        let dz = self.pos[2] - point[2];
        dx * dx + dy * dy + dz * dz
    }
}

/// An atom for SC calculation.
#[derive(Debug, Clone)]
struct ScAtom {
    coor: Vector3<f64>,
    radius: f64,
    molecule: usize,
    is_interface: bool,
}

/// Calculate shape complementarity between two chain groups.
pub fn get_sc(pdb: &PDB, groups: &str, model_num: usize, settings: Option<ScSettings>) -> ScResult {
    let settings = settings.unwrap_or_default();

    let pdb_prepared = prepare_pdb_for_sasa(pdb, true, true, "");
    let pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    let all_chains: HashSet<String> = pdb_filtered.chains().map(|c| c.id().to_string()).collect();
    let (group1_chains, group2_chains) = parse_groups(&all_chains, groups);

    let tree = pdb_filtered.create_hierarchy_rtree();
    let sep_sq = settings.separation * settings.separation;
    let rp = settings.probe_radius;

    let mut atoms: Vec<ScAtom> = Vec::new();
    let mut surface_results = [ScSurfaceResult::default(), ScSurfaceResult::default()];

    // Collect atoms
    for hier in pdb_filtered.atoms_with_hierarchy() {
        let chain_id = hier.chain().id().to_string();
        let molecule = if group1_chains.contains(&chain_id) {
            0
        } else if group2_chains.contains(&chain_id) {
            1
        } else {
            continue;
        };

        if hier.atom().element() == Some(&Element::H) {
            continue;
        }

        let vdw_radius = hier
            .atom()
            .element()
            .and_then(|e| e.atomic_radius().van_der_waals)
            .unwrap_or(1.5);

        let pos = hier.atom().pos();
        let coor = Vector3::new(pos.0, pos.1, pos.2);

        let other_group = if molecule == 0 { &group2_chains } else { &group1_chains };
        let is_interface = tree
            .locate_within_distance(pos, sep_sq)
            .any(|neighbor| other_group.contains(neighbor.chain().id()));

        surface_results[molecule].n_atoms += 1;
        if is_interface {
            surface_results[molecule].n_buried_atoms += 1;
        } else {
            surface_results[molecule].n_blocked_atoms += 1;
        }

        atoms.push(ScAtom { coor, radius: vdw_radius, molecule, is_interface });
    }

    if atoms.is_empty() || surface_results[0].n_atoms == 0 || surface_results[1].n_atoms == 0 {
        return ScResult {
            sc: 0.0, distance: 0.0, area: 0.0,
            surfaces: surface_results, valid: false,
        };
    }

    // Separate atoms by molecule
    let mol0_atoms: Vec<(Vector3<f64>, f64)> = atoms.iter()
        .filter(|a| a.molecule == 0)
        .map(|a| (a.coor, a.radius))
        .collect();
    let mol1_atoms: Vec<(Vector3<f64>, f64)> = atoms.iter()
        .filter(|a| a.molecule == 1)
        .map(|a| (a.coor, a.radius))
        .collect();

    // Generate surface dots (parallel)
    let dots_results: Vec<(usize, Vec<Dot>)> = atoms
        .par_iter()
        .enumerate()
        .filter_map(|(atom_idx, atom)| {
            if !atom.is_interface {
                return None;
            }

            let radius_i = atom.radius;
            let expanded_radius_i = atom.radius + rp;
            let surface_area = 4.0 * std::f64::consts::PI * radius_i * radius_i;
            let n_dots = (surface_area * settings.density).ceil() as usize;
            let area_per_dot = surface_area / n_dots as f64;

            let golden_ratio = (1.0 + 5.0_f64.sqrt()) / 2.0;
            let mut dots = Vec::new();

            let same_mol_atoms: Vec<(Vector3<f64>, f64)> = atoms.iter()
                .enumerate()
                .filter(|(i, a)| *i != atom_idx && a.molecule == atom.molecule)
                .map(|(_, a)| (a.coor, a.radius))
                .collect();

            let other_mol_atoms = if atom.molecule == 0 { &mol1_atoms } else { &mol0_atoms };

            for i in 0..n_dots {
                let theta = 2.0 * std::f64::consts::PI * (i as f64) / golden_ratio;
                let phi = (1.0 - 2.0 * (i as f64 + 0.5) / n_dots as f64).acos();

                let normal = Vector3::new(
                    phi.sin() * theta.cos(),
                    phi.sin() * theta.sin(),
                    phi.cos(),
                );

                let point = atom.coor + normal * radius_i;
                let pcen = atom.coor + normal * expanded_radius_i;

                // Collision check
                let mut collision = false;
                for (other_coor, other_radius) in &same_mol_atoms {
                    let other_expanded = other_radius + rp;
                    if (pcen - other_coor).norm_squared() < other_expanded * other_expanded {
                        collision = true;
                        break;
                    }
                }
                if collision {
                    continue;
                }

                // Burial check
                let mut buried = false;
                for (other_coor, other_radius) in other_mol_atoms {
                    let other_expanded = other_radius + rp;
                    if (pcen - other_coor).norm_squared() <= other_expanded * other_expanded {
                        buried = true;
                        break;
                    }
                }

                dots.push(Dot {
                    coor: point,
                    outnml: normal,
                    area: area_per_dot,
                    buried,
                });
            }

            if dots.is_empty() { None } else { Some((atom.molecule, dots)) }
        })
        .collect();

    let mut dots: [Vec<Dot>; 2] = [Vec::new(), Vec::new()];
    for (mol, atom_dots) in dots_results {
        dots[mol].extend(atom_dots);
    }

    surface_results[0].n_dots = dots[0].len();
    surface_results[1].n_dots = dots[1].len();

    if dots[0].is_empty() || dots[1].is_empty() {
        return ScResult {
            sc: 0.0, distance: 0.0, area: 0.0,
            surfaces: surface_results, valid: false,
        };
    }

    // Build RTree of accessible dots for fast trimming
    let band_sq = settings.band * settings.band;
    let mut trimmed_indices: [Vec<usize>; 2] = [Vec::new(), Vec::new()];

    for mol in 0..2 {
        let sdots = &dots[mol];
        
        // Build RTree of accessible (non-buried) dots
        let accessible_dots: Vec<DotPoint> = sdots.iter()
            .enumerate()
            .filter(|(_, d)| !d.buried)
            .map(|(idx, d)| DotPoint {
                idx,
                pos: [d.coor.x, d.coor.y, d.coor.z],
            })
            .collect();
        
        let accessible_tree = RTree::bulk_load(accessible_dots);
        
        // For each buried dot, check if any accessible dot is within band distance
        let indices: Vec<usize> = (0..sdots.len())
            .into_par_iter()
            .filter(|&idx| {
                let dot = &sdots[idx];
                if !dot.buried {
                    return false;
                }
                let pos = [dot.coor.x, dot.coor.y, dot.coor.z];
                // Use RTree to check for nearby accessible dots
                accessible_tree.locate_within_distance(pos, band_sq).next().is_none()
            })
            .collect();

        surface_results[mol].n_trimmed_dots = indices.len();
        surface_results[mol].trimmed_area = indices.iter().map(|&i| sdots[i].area).sum();
        trimmed_indices[mol] = indices;
    }

    if trimmed_indices[0].is_empty() || trimmed_indices[1].is_empty() {
        return ScResult {
            sc: 0.0, distance: 0.0, area: 0.0,
            surfaces: surface_results, valid: false,
        };
    }

    // Build RTree for neighbor finding
    let mut their_trees: [Option<RTree<DotPoint>>; 2] = [None, None];
    for mol in 0..2 {
        let their_dots: Vec<DotPoint> = trimmed_indices[mol].iter()
            .map(|&idx| DotPoint {
                idx,
                pos: [dots[mol][idx].coor.x, dots[mol][idx].coor.y, dots[mol][idx].coor.z],
            })
            .collect();
        their_trees[mol] = Some(RTree::bulk_load(their_dots));
    }

    // Calculate SC scores (parallel)
    for (my, their) in [(0, 1), (1, 0)] {
        let my_indices = &trimmed_indices[my];
        let my_dots = &dots[my];
        let their_dots = &dots[their];
        let their_tree = their_trees[their].as_ref().unwrap();

        let results: Vec<(f64, f64)> = my_indices
            .par_iter()
            .filter_map(|&idx| {
                let dot1 = &my_dots[idx];
                let pos = [dot1.coor.x, dot1.coor.y, dot1.coor.z];
                
                // Find nearest neighbor using RTree
                their_tree.nearest_neighbor(&pos).map(|nearest| {
                    let neighbor = &their_dots[nearest.idx];
                    let dist = (dot1.coor - neighbor.coor).norm();
                    let mut r = dot1.outnml.dot(&neighbor.outnml);
                    r *= (-dist * dist * settings.weight).exp();
                    r = r.clamp(-0.999, 0.999);
                    (dist, -r)
                })
            })
            .collect();

        if results.is_empty() {
            continue;
        }

        let mut distances: Vec<f64> = results.iter().map(|(d, _)| *d).collect();
        let mut scores: Vec<f64> = results.iter().map(|(_, s)| *s).collect();

        let n = distances.len() as f64;
        let d_mean = distances.iter().sum::<f64>() / n;
        let s_mean = -scores.iter().sum::<f64>() / n;

        let d_median = {
            let mid = distances.len() / 2;
            *distances.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap()).1
        };
        let s_median = {
            let mid = scores.len() / 2;
            *scores.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap()).1
        };

        surface_results[my].d_mean = d_mean;
        surface_results[my].d_median = d_median;
        surface_results[my].s_mean = s_mean;
        surface_results[my].s_median = s_median;
    }

    let sc = (surface_results[0].s_median + surface_results[1].s_median) / 2.0;
    let distance = (surface_results[0].d_median + surface_results[1].d_median) / 2.0;
    let area = surface_results[0].trimmed_area + surface_results[1].trimmed_area;

    ScResult {
        sc, distance, area,
        surfaces: surface_results,
        valid: true,
    }
}

/// Get the shape complementarity score as a simple f64 value.
pub fn get_sc_score(pdb: &PDB, groups: &str, model_num: usize) -> f64 {
    let result = get_sc(pdb, groups, model_num, None);
    if result.valid { result.sc } else { -1.0 }
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
    }

    #[test]
    fn test_sc_h_vs_l() {
        let pdb = load_multi_chain();
        let result = get_sc(&pdb, "H/L", 0, None);

        eprintln!("H vs L: sc={:.6}, expected~0.714", result.sc);
        eprintln!("  distance={:.4}, expected~0.56", result.distance);
        eprintln!("  Surface 0: n_dots={}, trimmed={}", result.surfaces[0].n_dots, result.surfaces[0].n_trimmed_dots);
        eprintln!("  Surface 1: n_dots={}, trimmed={}", result.surfaces[1].n_dots, result.surfaces[1].n_trimmed_dots);

        assert!(result.valid, "SC calculation should be valid");

        // sc-rs reference: 0.714. Allow ~0.15 tolerance for implementation differences
        let expected = 0.714;
        let tolerance = 0.15;
        assert!(
            (result.sc - expected).abs() < tolerance,
            "SC score {} should be within {} of expected {}",
            result.sc, tolerance, expected
        );
    }

    #[test]
    fn test_sc_h_vs_c() {
        let pdb = load_multi_chain();
        let result = get_sc(&pdb, "H/C", 0, None);

        eprintln!("H vs C: sc={:.6}, expected~0.785", result.sc);
        assert!(result.valid);

        // sc-rs reference: 0.785. Allow ~0.25 tolerance (larger interface, more variance)
        let expected = 0.785;
        let tolerance = 0.25;
        assert!((result.sc - expected).abs() < tolerance);
    }

    #[test]
    fn test_sc_l_vs_c() {
        let pdb = load_multi_chain();
        let result = get_sc(&pdb, "L/C", 0, None);

        eprintln!("L vs C: sc={:.6}, expected~0.734", result.sc);
        assert!(result.valid);

        // sc-rs reference: 0.734. Allow ~0.15 tolerance
        let expected = 0.734;
        let tolerance = 0.15;
        assert!((result.sc - expected).abs() < tolerance);
    }
}
