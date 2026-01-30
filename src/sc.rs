//! Shape Complementarity (SC) calculations.

use crate::sasa::{filter_pdb_by_model, prepare_pdb_for_sasa};
use crate::utils::parse_groups;
use nalgebra::Vector3;
use pdbtbx::*;
use rayon::prelude::*;
use std::collections::HashSet;

const GAUSSIAN_W: f64 = 0.5;
const PERIPHERAL_BAND: f64 = 1.5;
const DOT_DENSITY: f64 = 15.0;
const PROBE_RADIUS: f64 = 1.7;
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

#[derive(Debug, Clone)]
struct Dot {
    coor: Vector3<f64>,
    outnml: Vector3<f64>,
    area: f64,
    buried: bool,
}

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

    let mut atoms: Vec<ScAtom> = Vec::new();
    let mut surface_results = [ScSurfaceResult::default(), ScSurfaceResult::default()];

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

    let rp = settings.probe_radius;
    
    let mol0_atoms: Vec<&ScAtom> = atoms.iter().filter(|a| a.molecule == 0).collect();
    let mol1_atoms: Vec<&ScAtom> = atoms.iter().filter(|a| a.molecule == 1).collect();

    let mut dots: [Vec<Dot>; 2] = [Vec::new(), Vec::new()];
    
    for (atom_idx, atom) in atoms.iter().enumerate() {
        if !atom.is_interface {
            continue;
        }
        
        let expanded_radius = atom.radius + rp;
        let surface_area = 4.0 * std::f64::consts::PI * expanded_radius * expanded_radius;
        let n_dots = (surface_area * settings.density).ceil() as usize;
        let area_per_dot = surface_area / n_dots as f64;

        let golden_ratio = (1.0 + 5.0_f64.sqrt()) / 2.0;
        
        for i in 0..n_dots {
            let theta = 2.0 * std::f64::consts::PI * (i as f64) / golden_ratio;
            let phi = (1.0 - 2.0 * (i as f64 + 0.5) / n_dots as f64).acos();
            
            let outnml = Vector3::new(
                phi.sin() * theta.cos(),
                phi.sin() * theta.sin(),
                phi.cos(),
            );
            let coor = atom.coor + outnml * expanded_radius;

            // Check collision with same-molecule atoms
            let mut collision = false;
            for (other_idx, other) in atoms.iter().enumerate() {
                if other_idx == atom_idx || other.molecule != atom.molecule {
                    continue;
                }
                let other_expanded = other.radius + rp;
                if (coor - other.coor).norm_squared() < other_expanded * other_expanded {
                    collision = true;
                    break;
                }
            }
            if collision {
                continue;
            }

            // Check burial by opposite molecule
            let other_mol_atoms = if atom.molecule == 0 { &mol1_atoms } else { &mol0_atoms };
            let mut buried = false;
            for other in other_mol_atoms.iter() {
                let other_expanded = other.radius + rp;
                if (coor - other.coor).norm_squared() <= other_expanded * other_expanded {
                    buried = true;
                    break;
                }
            }

            dots[atom.molecule].push(Dot { coor, outnml, area: area_per_dot, buried });
        }
    }

    surface_results[0].n_dots = dots[0].len();
    surface_results[1].n_dots = dots[1].len();

    if dots[0].is_empty() || dots[1].is_empty() {
        return ScResult {
            sc: 0.0, distance: 0.0, area: 0.0,
            surfaces: surface_results, valid: false,
        };
    }

    // Trim peripheral band - keep buried dots that are far from accessible dots
    // (i.e., the stable interior of the interface)
    let band_sq = settings.band * settings.band;
    let mut trimmed_indices: [Vec<usize>; 2] = [Vec::new(), Vec::new()];

    for mol in 0..2 {
        let sdots = &dots[mol];
        for (idx, dot) in sdots.iter().enumerate() {
            if !dot.buried {
                continue;
            }
            // Keep if no accessible dot within band distance (interior dots)
            let has_nearby_accessible = sdots.iter().enumerate().any(|(j, dot2)| {
                j != idx && !dot2.buried && (dot.coor - dot2.coor).norm_squared() <= band_sq
            });
            if !has_nearby_accessible {
                trimmed_indices[mol].push(idx);
            }
        }
        surface_results[mol].n_trimmed_dots = trimmed_indices[mol].len();
        surface_results[mol].trimmed_area = trimmed_indices[mol]
            .iter()
            .map(|&i| dots[mol][i].area)
            .sum();
    }

    if trimmed_indices[0].is_empty() || trimmed_indices[1].is_empty() {
        return ScResult {
            sc: 0.0, distance: 0.0, area: 0.0,
            surfaces: surface_results, valid: false,
        };
    }

    // Calculate SC scores using ALL buried dots (not just trimmed) for neighbor finding
    // but only calculate stats for trimmed dots
    for (my, their) in [(0, 1), (1, 0)] {
        let my_indices = &trimmed_indices[my];
        let my_dots = &dots[my];
        let their_dots = &dots[their];
        
        // Use ALL buried dots from their surface for neighbor search
        let their_buried: Vec<usize> = (0..their_dots.len())
            .filter(|&i| their_dots[i].buried)
            .collect();

        let results: Vec<(f64, f64)> = my_indices
            .par_iter()
            .filter_map(|&idx| {
                let dot1 = &my_dots[idx];
                
                let mut min_dist_sq = f64::MAX;
                let mut nearest: Option<&Dot> = None;
                
                // Find nearest BURIED neighbor from their surface
                for &idx2 in &their_buried {
                    let dot2 = &their_dots[idx2];
                    let d_sq = (dot1.coor - dot2.coor).norm_squared();
                    if d_sq < min_dist_sq {
                        min_dist_sq = d_sq;
                        nearest = Some(dot2);
                    }
                }

                nearest.map(|neighbor| {
                    let dist = min_dist_sq.sqrt();
                    let r = dot1.outnml.dot(&neighbor.outnml);
                    let weighted = r * (-dist * dist * settings.weight).exp();
                    let clamped = weighted.clamp(-0.999, 0.999);
                    (dist, -clamped)
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
        let s_mean = scores.iter().sum::<f64>() / n;

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
        
        eprintln!("H vs L: sc={:.6}, expected=0.725945", result.sc);
        eprintln!("  Surface 0: atoms={}, interface={}, dots={}, trimmed={}", 
            result.surfaces[0].n_atoms, result.surfaces[0].n_buried_atoms,
            result.surfaces[0].n_dots, result.surfaces[0].n_trimmed_dots);
        eprintln!("  Surface 1: atoms={}, interface={}, dots={}, trimmed={}", 
            result.surfaces[1].n_atoms, result.surfaces[1].n_buried_atoms,
            result.surfaces[1].n_dots, result.surfaces[1].n_trimmed_dots);
        eprintln!("  s_median: [{:.4}, {:.4}]", 
            result.surfaces[0].s_median, result.surfaces[1].s_median);
        eprintln!("  d_median: [{:.4}, {:.4}]", 
            result.surfaces[0].d_median, result.surfaces[1].d_median);
        
        assert!(result.valid, "SC calculation should be valid");
        
        let expected = 0.725945;
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
        
        eprintln!("H vs C: sc={:.6}, expected=0.778542", result.sc);
        assert!(result.valid);
        
        let expected = 0.778542;
        let tolerance = 0.15;
        assert!((result.sc - expected).abs() < tolerance);
    }

    #[test]
    fn test_sc_l_vs_c() {
        let pdb = load_multi_chain();
        let result = get_sc(&pdb, "L/C", 0, None);
        
        eprintln!("L vs C: sc={:.6}, expected=0.733599", result.sc);
        assert!(result.valid);
        
        let expected = 0.733599;
        let tolerance = 0.15;
        assert!((result.sc - expected).abs() < tolerance);
    }
}
