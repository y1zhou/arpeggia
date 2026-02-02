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

/// A surface dot with position and normal vector.
#[derive(Debug, Clone)]
struct Dot {
    coor: Vector3<f64>,
    outnml: Vector3<f64>,
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
/// Returns the SC score as f64, or -1.0 if calculation fails.
pub fn get_sc(pdb: &PDB, groups: &str, model_num: usize) -> f64 {
    let rp = PROBE_RADIUS;
    let density = DOT_DENSITY;
    let band = PERIPHERAL_BAND;
    let sep = SEPARATION_CUTOFF;
    let weight = GAUSSIAN_W;

    let pdb_prepared = prepare_pdb_for_sasa(pdb, true, true, "");
    let mut pdb_filtered = filter_pdb_by_model(&pdb_prepared, model_num);

    let all_chains: HashSet<String> = pdb_filtered.chains().map(|c| c.id().to_string()).collect();
    let (group1_chains, group2_chains) = parse_groups(&all_chains, groups);

    // Remove chains not in either group
    pdb_filtered.remove_chains_by(|c| {
        let cid = c.id().to_string();
        !(group1_chains.contains(&cid) || group2_chains.contains(&cid))
    });

    let tree = pdb_filtered.create_hierarchy_rtree();
    let sep_sq = sep * sep;

    // Collect atoms using filter_map/collect pattern
    let atoms: Vec<ScAtom> = pdb_filtered
        .atoms_with_hierarchy()
        .filter_map(|hier| {
            let chain_id = hier.chain().id().to_string();
            let molecule = if group1_chains.contains(&chain_id) {
                0
            } else if group2_chains.contains(&chain_id) {
                1
            } else {
                return None;
            };

            // Skip hydrogens
            if hier.atom().element() == Some(&Element::H) {
                return None;
            }

            let vdw_radius = hier
                .atom()
                .element()
                .and_then(|e| e.atomic_radius().van_der_waals)
                .unwrap_or(1.5);

            let pos = hier.atom().pos();
            let coor = Vector3::new(pos.0, pos.1, pos.2);

            // Check if interface atom (near opposite molecule)
            let other_group = if molecule == 0 { &group2_chains } else { &group1_chains };
            let is_interface = tree
                .locate_within_distance(pos, sep_sq)
                .any(|neighbor| other_group.contains(neighbor.chain().id()));

            Some(ScAtom {
                coor,
                radius: vdw_radius,
                molecule,
                is_interface,
            })
        })
        .collect();

    if atoms.is_empty() {
        return -1.0;
    }

    let mol0_count = atoms.iter().filter(|a| a.molecule == 0).count();
    let mol1_count = atoms.iter().filter(|a| a.molecule == 1).count();
    if mol0_count == 0 || mol1_count == 0 {
        return -1.0;
    }

    // Separate atoms by molecule for efficient checks
    let mol0_atoms: Vec<(Vector3<f64>, f64)> = atoms.iter()
        .filter(|a| a.molecule == 0)
        .map(|a| (a.coor, a.radius))
        .collect();
    let mol1_atoms: Vec<(Vector3<f64>, f64)> = atoms.iter()
        .filter(|a| a.molecule == 1)
        .map(|a| (a.coor, a.radius))
        .collect();

    // Generate surface dots using Fibonacci sphere sampling (parallel)
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
            let n_dots = (surface_area * density).ceil() as usize;

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

                // Collision check with same-molecule atoms
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

                // Burial check against opposite molecule
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

    if dots[0].is_empty() || dots[1].is_empty() {
        return -1.0;
    }

    // Trim peripheral band using RTree for efficiency
    let band_sq = band * band;
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
        
        // Keep buried dots that have no accessible dots within band distance
        let indices: Vec<usize> = (0..sdots.len())
            .into_par_iter()
            .filter(|&idx| {
                let dot = &sdots[idx];
                if !dot.buried {
                    return false;
                }
                let pos = [dot.coor.x, dot.coor.y, dot.coor.z];
                accessible_tree.locate_within_distance(pos, band_sq).next().is_none()
            })
            .collect();

        trimmed_indices[mol] = indices;
    }

    if trimmed_indices[0].is_empty() || trimmed_indices[1].is_empty() {
        return -1.0;
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

    // Calculate SC scores
    let mut s_medians = [0.0f64; 2];

    for (my, their) in [(0, 1), (1, 0)] {
        let my_indices = &trimmed_indices[my];
        let my_dots = &dots[my];
        let their_dots = &dots[their];
        let their_tree = their_trees[their].as_ref().unwrap();

        let scores: Vec<f64> = my_indices
            .par_iter()
            .filter_map(|&idx| {
                let dot1 = &my_dots[idx];
                let pos = [dot1.coor.x, dot1.coor.y, dot1.coor.z];
                
                // Find nearest neighbor using RTree
                their_tree.nearest_neighbor(&pos).map(|nearest| {
                    let neighbor = &their_dots[nearest.idx];
                    let dist = (dot1.coor - neighbor.coor).norm();
                    let mut r = dot1.outnml.dot(&neighbor.outnml);
                    r *= (-dist * dist * weight).exp();
                    r.clamp(-0.999, 0.999)
                })
            })
            .collect();

        if scores.is_empty() {
            return -1.0;
        }

        // Calculate median
        // The dot product of opposite-facing (complementary) normals is negative.
        // We negate to get positive SC values for well-fitting surfaces.
        let mut scores = scores;
        let mid = scores.len() / 2;
        let s_median = *scores.select_nth_unstable_by(mid, |a, b| a.partial_cmp(b).unwrap()).1;
        s_medians[my] = -s_median;
    }

    (s_medians[0] + s_medians[1]) / 2.0
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
    fn test_sc_h_vs_l() {
        let pdb = load_multi_chain();
        let sc = get_sc(&pdb, "H/L", 0);

        eprintln!("H vs L: sc={:.6}, expected~0.714 (sc-rs reference)", sc);

        assert!(sc > 0.0, "SC calculation should succeed");

        // Target: within ~0.15 of sc-rs reference (0.714)
        let expected = 0.714;
        let tolerance = 0.15;
        assert!(
            (sc - expected).abs() < tolerance,
            "SC score {} should be within {} of expected {}",
            sc, tolerance, expected
        );
    }

    #[test]
    fn test_sc_h_vs_c() {
        let pdb = load_multi_chain();
        let sc = get_sc(&pdb, "H/C", 0);

        eprintln!("H vs C: sc={:.6}, expected~0.785 (sc-rs reference)", sc);
        assert!(sc > 0.0);

        let expected = 0.785;
        let tolerance = 0.25;
        assert!((sc - expected).abs() < tolerance);
    }

    #[test]
    fn test_sc_l_vs_c() {
        let pdb = load_multi_chain();
        let sc = get_sc(&pdb, "L/C", 0);

        eprintln!("L vs C: sc={:.6}, expected~0.734 (sc-rs reference)", sc);
        assert!(sc > 0.0);

        let expected = 0.734;
        let tolerance = 0.15;
        assert!((sc - expected).abs() < tolerance);
    }
}
