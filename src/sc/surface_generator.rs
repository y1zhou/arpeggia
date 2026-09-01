//! Connolly molecular surface generator for SC calculations.
//! Ported from <https://github.com/cytokineking/sc-rs>

use core::f64;
use std::f64::consts::PI;
use std::fmt;

use super::atomic_radii::{EMBEDDED_ATOMIC_RADII, wildcard_match};
use super::types::*;
use super::vector3::Vec3;
use super::{DOT_DENSITY, PROBE_RADIUS};
use rayon::prelude::*;

/// Error type for surface calculation operations.
#[derive(Debug)]
pub(super) enum SurfaceCalculatorError {
    /// No atoms defined
    NoAtoms,
    /// A selected chain group contains no usable atoms.
    NoAtomsForGroup(usize),
    /// An atom has no usable radius.
    MissingRadius(String),
    /// Overlapping atoms detected
    Coincident(String),
    /// The selected groups do not form a sampled interface.
    NoInterface,
    /// Invalid public calculation arguments.
    InvalidInput(String),
}

impl fmt::Display for SurfaceCalculatorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SurfaceCalculatorError::NoAtoms => write!(f, "No atoms defined"),
            SurfaceCalculatorError::NoAtomsForGroup(group) => {
                write!(f, "No atoms for chain group {group}")
            }
            SurfaceCalculatorError::MissingRadius(message) => f.write_str(message),
            SurfaceCalculatorError::Coincident(msg) => {
                write!(f, "Overlapping atoms detected: {msg}")
            }
            SurfaceCalculatorError::NoInterface => {
                write!(f, "selected groups do not form a sampled interface")
            }
            SurfaceCalculatorError::InvalidInput(message) => write!(f, "Invalid input: {message}"),
        }
    }
}

impl std::error::Error for SurfaceCalculatorError {}

#[derive(Default)]
pub(super) struct SurfaceGenerator {
    pub(super) run: RunState,
}

#[derive(Clone, Default)]
pub(super) struct RunState {
    pub(super) atoms: Vec<ScAtom>,
    pub(super) probes: Vec<Probe>,
    pub(super) dots: [Vec<Dot>; 2],
    pub(super) trimmed_dots: [Vec<usize>; 2],
    pub(super) results: Results,
}

impl SurfaceGenerator {
    pub(super) fn get_atom_radius(resn: &str, atomn: &str) -> Result<f64, SurfaceCalculatorError> {
        let resn = if crate::contacts::is_histidine(resn) {
            "HIS"
        } else {
            resn
        };
        // First try the embedded radii table
        for radius in EMBEDDED_ATOMIC_RADII {
            if !wildcard_match(resn, radius.residue) {
                continue;
            }
            if !wildcard_match(atomn, radius.atom) {
                continue;
            }
            return Ok(radius.radius);
        }
        Err(SurfaceCalculatorError::MissingRadius(format!(
            "no radius for {resn}:{atomn}"
        )))
    }

    pub(super) fn generate_molecular_surfaces(&mut self) -> Result<(), SurfaceCalculatorError> {
        if self.run.atoms.is_empty() {
            return Err(SurfaceCalculatorError::NoAtoms);
        }
        let len = self.run.atoms.len();
        for i in 0..len {
            let att = self.run.atoms[i].attention;
            if matches!(att, Attention::Far) {
                continue;
            }
            self.build_probes(i)?;
        }
        self.generate_contact_surface()?;
        self.generate_concave_surface()?;
        Ok(())
    }

    fn generate_contact_surface(&mut self) -> Result<(), SurfaceCalculatorError> {
        let rp = PROBE_RADIUS;
        let atoms: &Vec<ScAtom> = &self.run.atoms;
        let results: Vec<(usize, Vec<Dot>)> = (0..atoms.len())
            .into_par_iter()
            .filter_map(|i| {
                let a_i = &atoms[i];
                let att = a_i.attention;
                if matches!(att, Attention::Far) {
                    return None;
                }
                if !a_i.accessible {
                    return None;
                }
                let neighbors = &a_i.neighbor_indices;
                let mut north_dir = Vec3::new(0.0, 0.0, 1.0);
                let mut south_dir = Vec3::new(0.0, 0.0, -1.0);
                let mut equatorial_vector = Vec3::new(1.0, 0.0, 0.0);
                let radius_i = a_i.radius;
                let expanded_radius_i = a_i.radius + rp;
                if !neighbors.is_empty() {
                    let neighbor = &atoms[neighbors[0]];
                    north_dir = a_i.coor - neighbor.coor;
                    north_dir.normalize();
                    let mut temp_vec = Vec3::new(
                        north_dir.y * north_dir.y + north_dir.z * north_dir.z,
                        north_dir.x * north_dir.x + north_dir.z * north_dir.z,
                        north_dir.x * north_dir.x + north_dir.y * north_dir.y,
                    );
                    temp_vec.normalize();
                    let dt = temp_vec.dot(north_dir);
                    if dt.abs() > 0.99 {
                        temp_vec = Vec3::new(1.0, 0.0, 0.0);
                    }
                    equatorial_vector = north_dir.cross(temp_vec);
                    equatorial_vector.normalize();
                    let radius_neighbor = neighbor.radius;
                    let expanded_radius_j = neighbor.radius + rp;
                    let dij = a_i.coor.distance(neighbor.coor);
                    let unit_axis = (neighbor.coor - a_i.coor) / dij;
                    let asymmetry_term = (expanded_radius_i * expanded_radius_i
                        - expanded_radius_j * expanded_radius_j)
                        / dij;
                    let midplane_center =
                        (a_i.coor + neighbor.coor) * 0.5 + (unit_axis * (asymmetry_term * 0.5));
                    let mut far_term = (expanded_radius_i + expanded_radius_j)
                        * (expanded_radius_i + expanded_radius_j)
                        - dij * dij;
                    if far_term <= 0.0 {
                        return None;
                    }
                    far_term = far_term.sqrt();
                    let mut contain_term = dij * dij - (radius_i - radius_neighbor).powi(2);
                    if contain_term <= 0.0 {
                        return None;
                    }
                    contain_term = contain_term.sqrt();
                    let ring_radius = 0.5 * far_term * contain_term / dij;
                    let ring_point =
                        midplane_center + (equatorial_vector.cross(north_dir) * ring_radius);
                    south_dir = (ring_point - a_i.coor) / expanded_radius_i;
                    if north_dir.cross(south_dir).dot(equatorial_vector) <= 0.0 {
                        return None;
                    }
                }
                let mut lats: Vec<Vec3> = Vec::new();
                let o = Vec3::zero();
                geom_sample_arc(
                    o,
                    radius_i,
                    equatorial_vector,
                    DOT_DENSITY,
                    north_dir,
                    south_dir,
                    &mut lats,
                )
                .ok()?;
                if lats.is_empty() {
                    return None;
                }
                let mut dots: Vec<Dot> = Vec::new();
                let mut points: Vec<Vec3> = Vec::new();
                for ilat in &lats {
                    let dt = ilat.dot(north_dir);
                    let cen = a_i.coor + (north_dir * dt);
                    let mut rad = radius_i * radius_i - dt * dt;
                    if rad <= 0.0 {
                        continue;
                    }
                    rad = rad.sqrt();
                    points.clear();
                    geom_sample_circle(cen, rad, north_dir, DOT_DENSITY, &mut points).ok()?;
                    if points.is_empty() {
                        continue;
                    }
                    for &point in &points {
                        let pcen = a_i.coor + ((point - a_i.coor) * (expanded_radius_i / radius_i));
                        // collision with same-molecule neighbors (skip first neighbor)
                        let mut coll = false;
                        for &idx in neighbors.iter().skip(1) {
                            let a = &atoms[idx];
                            if pcen.distance(a.coor) <= (a.radius + rp) {
                                coll = true;
                                break;
                            }
                        }
                        if coll {
                            continue;
                        }
                        // burial check against opposite molecule
                        let other_mol = u8::from(a_i.molecule == 0);
                        let mut buried = false;
                        for b in atoms {
                            if b.molecule != other_mol {
                                continue;
                            }
                            let erl = b.radius + rp;
                            let d = pcen.distance_squared(b.coor);
                            if d <= erl * erl {
                                buried = true;
                                break;
                            }
                        }
                        let outnml = if rp <= 0.0 {
                            point - a_i.coor
                        } else {
                            (pcen - point) / rp
                        };
                        dots.push(Dot {
                            coor: point,
                            outnml,
                            buried,
                        });
                    }
                }
                if dots.is_empty() {
                    None
                } else {
                    Some((usize::from(a_i.molecule), dots))
                }
            })
            .collect();
        for (mol, mut dots) in results {
            self.run.dots[mol].append(&mut dots);
        }
        Ok(())
    }

    fn build_probes(&mut self, atom_index: usize) -> Result<(), SurfaceCalculatorError> {
        let atom1_coor: Vec3;
        let expanded_radius_i: f64;
        let neighbor_indices: Vec<usize>;
        {
            let atom1 = &self.run.atoms[atom_index];
            atom1_coor = atom1.coor;
            expanded_radius_i = atom1.radius + PROBE_RADIUS;
            neighbor_indices = atom1.neighbor_indices.clone();
        }
        let num_neighbors = neighbor_indices.len();

        for &j in &neighbor_indices {
            let atom2 = &self.run.atoms[j];
            if atom2.atomi <= self.run.atoms[atom_index].atomi {
                continue;
            }
            let expanded_radius_j = atom2.radius + PROBE_RADIUS;
            let ij_dist2 = self.run.atoms[atom_index]
                .neighbor_distance_squared(j)
                .unwrap();
            let dist_ij = ij_dist2.sqrt();

            let unit_axis = (atom2.coor - atom1_coor) / dist_ij;
            let asymmetry_term = (expanded_radius_i * expanded_radius_i
                - expanded_radius_j * expanded_radius_j)
                / dist_ij;
            let midplane_center =
                (atom1_coor + atom2.coor) * 0.5 + (unit_axis * (asymmetry_term * 0.5));
            let mut far_term = (expanded_radius_i + expanded_radius_j)
                * (expanded_radius_i + expanded_radius_j)
                - ij_dist2;
            if far_term <= 0.0 {
                continue;
            }
            far_term = far_term.sqrt();
            let mut contain_term =
                ij_dist2 - (self.run.atoms[atom_index].radius - atom2.radius).powi(2);
            if contain_term <= 0.0 {
                continue;
            }
            contain_term = contain_term.sqrt();
            let ring_radius = 0.5 * far_term * contain_term / dist_ij;
            if num_neighbors <= 1 {
                self.run.atoms[atom_index].accessible = true;
                self.run.atoms[j].accessible = true;
                break;
            }
            self.build_probe_triplets(atom_index, j, unit_axis, midplane_center, ring_radius)?;
            let has_point_cusp = asymmetry_term.abs() < dist_ij;
            let atom2_att = self.run.atoms[j].attention;
            if !matches!(self.run.atoms[atom_index].attention, Attention::Far)
                || !matches!(atom2_att, Attention::Far)
            {
                self.emit_reentrant_surface(
                    atom_index,
                    j,
                    unit_axis,
                    midplane_center,
                    ring_radius,
                    has_point_cusp,
                )?;
            }
        }
        Ok(())
    }

    fn build_probe_triplets(
        &mut self,
        atom1_index: usize,
        atom2_index: usize,
        unit_axis: Vec3,
        midplane_center: Vec3,
        ring_radius: f64,
    ) -> Result<(), SurfaceCalculatorError> {
        let atom1_coor = self.run.atoms[atom1_index].coor;
        let neighbor_indices = self.run.atoms[atom1_index].neighbor_indices.clone();
        let expanded_radius_i = self.run.atoms[atom1_index].radius + PROBE_RADIUS;

        let atom2 = &self.run.atoms[atom2_index];
        let expanded_radius_j = atom2.radius + PROBE_RADIUS;
        let atom2_natom = atom2.atomi;
        let atom2_att = atom2.attention;
        let mut made_probe = false;
        for &k in &neighbor_indices {
            let atom3 = &self.run.atoms[k];
            if atom3.atomi <= atom2_natom {
                continue;
            }
            let expanded_radius_k = atom3.radius + PROBE_RADIUS;

            // atom1 neighbor atom3 is not necessarily atom2 neighbor
            match atom2.neighbor_distance_squared(k) {
                None => continue,
                Some(dist_jk2) => {
                    if dist_jk2.sqrt() >= expanded_radius_j + expanded_radius_k {
                        continue;
                    }
                }
            }
            let dist_ik = self.run.atoms[atom1_index]
                .neighbor_distance_squared(k)
                .unwrap()
                .sqrt();
            if dist_ik >= expanded_radius_i + expanded_radius_k {
                continue;
            }
            if matches!(self.run.atoms[atom1_index].attention, Attention::Far)
                && matches!(atom2_att, Attention::Far)
                && matches!(atom3.attention, Attention::Far)
            {
                continue;
            }
            let unit_axis_ik = (atom3.coor - atom1_coor) / dist_ik;
            let wedge_angle = unit_axis.dot(unit_axis_ik).acos();
            let sin_wedge = wedge_angle.sin();
            if sin_wedge <= 0.0 {
                let dtijk2 = midplane_center.distance(atom3.coor);
                let rkp2 = expanded_radius_k * expanded_radius_k - ring_radius * ring_radius;
                if dtijk2 < rkp2 {
                    return Ok(());
                }
                continue;
            }
            let axis_normal = unit_axis.cross(unit_axis_ik) / sin_wedge;
            let perp_tangent = axis_normal.cross(unit_axis);
            let asymmetry_term_ik = (expanded_radius_i * expanded_radius_i
                - expanded_radius_k * expanded_radius_k)
                / dist_ik;
            let midpoint_ik =
                (atom1_coor + atom3.coor) * 0.5 + unit_axis_ik * (asymmetry_term_ik * 0.5);
            let mut componentwise = midpoint_ik - midplane_center;
            componentwise = Vec3::new(
                unit_axis_ik.x * componentwise.x,
                unit_axis_ik.y * componentwise.y,
                unit_axis_ik.z * componentwise.z,
            );
            let component_sum = componentwise.x + componentwise.y + componentwise.z;
            let torus_center = midplane_center + perp_tangent * (component_sum / sin_wedge);
            let mut height =
                expanded_radius_i * expanded_radius_i - torus_center.distance_squared(atom1_coor);
            if height <= 0.0 {
                continue;
            }
            height = height.sqrt();
            for is0 in 1..=2 {
                let sign_choice = 3 - 2 * is0;
                let probe_center = torus_center + axis_normal * (height * f64::from(sign_choice));
                if self.check_atom_collision2_idx(probe_center, atom2_index, k, &neighbor_indices) {
                    continue;
                }
                let mut probe = Probe {
                    atom_indices: [0; 3],
                    height,
                    point: probe_center,
                    alt: axis_normal * f64::from(sign_choice),
                };
                if sign_choice > 0 {
                    probe.atom_indices = [atom1_index, atom2_index, k];
                } else {
                    probe.atom_indices = [atom2_index, atom1_index, k];
                }
                self.run.probes.push(probe);
                made_probe = true;
            }
        }
        if made_probe {
            self.run.atoms[atom1_index].accessible = true;
        }
        Ok(())
    }

    fn emit_reentrant_surface(
        &mut self,
        atom1_index: usize,
        atom2_index: usize,
        unit_axis: Vec3,
        midplane_center: Vec3,
        ring_radius: f64,
        has_point_cusp: bool,
    ) -> Result<(), SurfaceCalculatorError> {
        let neighbors = self.run.atoms[atom1_index].neighbor_indices.clone();
        let density = DOT_DENSITY;
        let expanded_radius_i = self.run.atoms[atom1_index].radius + PROBE_RADIUS;
        let expanded_radius_j = self.run.atoms[atom2_index].radius + PROBE_RADIUS;
        let roll_circle_radius_i =
            ring_radius * self.run.atoms[atom1_index].radius / expanded_radius_i;
        let roll_circle_radius_j =
            ring_radius * self.run.atoms[atom2_index].radius / expanded_radius_j;
        let mut belt_radius = ring_radius - PROBE_RADIUS;
        if belt_radius <= 0.0 {
            belt_radius = 0.0;
        }
        let mean_radius = (roll_circle_radius_i + 2.0 * belt_radius + roll_circle_radius_j) / 4.0;
        let eccentricity = mean_radius / ring_radius;
        let effective_density = eccentricity * eccentricity * density;
        let mut subs: Vec<Vec3> = Vec::new();
        geom_sample_circle(
            midplane_center,
            ring_radius,
            unit_axis,
            effective_density,
            &mut subs,
        )?;
        if subs.is_empty() {
            return Ok(());
        }
        let atom2_atomi = self.run.atoms[atom2_index].atomi;
        for ring_point in subs {
            let mut tooclose = false;
            for &ni in &neighbors {
                let neighbor = &self.run.atoms[ni];
                if neighbor.atomi == atom2_atomi {
                    continue;
                }
                let expanded_neighbor_radius = (neighbor.radius + PROBE_RADIUS).powi(2);
                let d2 = ring_point.distance_squared(neighbor.coor);
                if d2 < expanded_neighbor_radius {
                    tooclose = true;
                    break;
                }
            }
            if tooclose {
                continue;
            }
            self.run.atoms[atom1_index].accessible = true;
            self.run.atoms[atom2_index].accessible = true;
            let vec_pi = (self.run.atoms[atom1_index].coor - ring_point) / expanded_radius_i;
            let vec_pj = (self.run.atoms[atom2_index].coor - ring_point) / expanded_radius_j;
            let mut toroid_axis = vec_pi.cross(vec_pj);
            toroid_axis.normalize();
            let mut cusp_term = PROBE_RADIUS * PROBE_RADIUS - ring_radius * ring_radius;
            let has_cusp_point = cusp_term > 0.0 && has_point_cusp;
            let (arc_end_i, arc_end_j) = if has_cusp_point {
                cusp_term = cusp_term.sqrt();
                let qij = midplane_center - unit_axis * cusp_term;
                (((qij - ring_point) / PROBE_RADIUS), Vec3::zero())
            } else {
                let mut pq = vec_pi + vec_pj;
                pq.normalize();
                (pq, pq)
            };
            let mut dot_tmp = arc_end_i.dot(vec_pi);
            if dot_tmp >= 1.0 || dot_tmp <= -1.0 {
                return Ok(());
            }
            dot_tmp = arc_end_j.dot(vec_pj);
            if dot_tmp >= 1.0 || dot_tmp <= -1.0 {
                return Ok(());
            }
            if !matches!(self.run.atoms[atom1_index].attention, Attention::Far) {
                let mut points: Vec<Vec3> = Vec::new();
                geom_sample_arc(
                    ring_point,
                    PROBE_RADIUS,
                    toroid_axis,
                    density,
                    vec_pi,
                    arc_end_i,
                    &mut points,
                )?;
                for &point in &points {
                    let molecule = self.run.atoms[atom1_index].molecule;
                    self.add_dot(molecule, point, ring_point);
                }
            }
            if matches!(self.run.atoms[atom2_index].attention, Attention::Far) {
                continue;
            }
            let mut points: Vec<Vec3> = Vec::new();
            geom_sample_arc(
                ring_point,
                PROBE_RADIUS,
                toroid_axis,
                density,
                arc_end_j,
                vec_pj,
                &mut points,
            )?;
            for point in points {
                let molecule2 = self.run.atoms[atom2_index].molecule;
                self.add_dot(molecule2, point, ring_point);
            }
        }
        Ok(())
    }

    fn check_atom_collision2_idx(
        &self,
        probe_center: Vec3,
        atom1_index: usize,
        atom2_index: usize,
        neighbor_indices: &[usize],
    ) -> bool {
        let atom1_natom = self.run.atoms[atom1_index].atomi;
        let atom2_natom = self.run.atoms[atom2_index].atomi;
        for &ni in neighbor_indices {
            let neighbor = &self.run.atoms[ni];
            if neighbor.atomi == atom1_natom || neighbor.atomi == atom2_natom {
                continue;
            }
            if probe_center.distance_squared(neighbor.coor)
                <= (neighbor.radius + PROBE_RADIUS).powi(2)
            {
                return true;
            }
        }
        false
    }

    fn generate_concave_surface(&mut self) -> Result<(), SurfaceCalculatorError> {
        let rp = PROBE_RADIUS;
        let rp2 = rp * rp;
        let atoms: &Vec<ScAtom> = &self.run.atoms;
        let probes: &Vec<Probe> = &self.run.probes;
        if probes.is_empty() {
            return Ok(());
        }
        let mut lowprobs: Vec<usize> = Vec::new();
        for (idx, probe) in probes.iter().enumerate() {
            if probe.height < rp {
                lowprobs.push(idx);
            }
        }
        let results: Vec<(Vec<Dot>, Vec<Dot>)> = (0..probes.len())
            .into_par_iter()
            .filter_map(|i| {
                let probe = &probes[i];
                let aidx = probe.atom_indices;
                let pijk = probe.point;
                let uijk = probe.alt;
                let hijk = probe.height;
                let density = DOT_DENSITY;
                let mut nears: Vec<usize> = Vec::new();
                for &lp in &lowprobs {
                    if lp == i {
                        continue;
                    }
                    let d2 = pijk.distance_squared(probes[lp].point);
                    if d2 <= 4.0 * rp2 {
                        nears.push(lp);
                    }
                }
                let mut vp = [Vec3::zero(); 3];
                for k in 0..3 {
                    vp[k] = atoms[aidx[k]].coor - pijk;
                    vp[k].normalize();
                }
                let mut vectors = [Vec3::zero(); 3];
                vectors[0] = vp[0].cross(vp[1]).normalized();
                vectors[1] = vp[1].cross(vp[2]).normalized();
                vectors[2] = vp[2].cross(vp[0]).normalized();
                let mut dm = -1.0;
                let mut mm = 0usize;
                for (k, vp_k) in vp.iter().enumerate() {
                    let dt = uijk.dot(*vp_k);
                    if dt > dm {
                        dm = dt;
                        mm = k;
                    }
                }
                let south_dir = uijk * -1.0;
                let mut arc_axis = vp[mm].cross(south_dir);
                arc_axis.normalize();
                let mut lats: Vec<Vec3> = Vec::new();
                let o = Vec3::zero();
                geom_sample_arc(o, rp, arc_axis, density, vp[mm], south_dir, &mut lats).ok()?;
                if lats.is_empty() {
                    return None;
                }
                let mut d0: Vec<Dot> = Vec::new();
                let mut d1: Vec<Dot> = Vec::new();
                let mut points: Vec<Vec3> = Vec::new();
                for ilat in &lats {
                    let dt = ilat.dot(south_dir);
                    let cen = south_dir * dt;
                    let mut rad = rp2 - dt * dt;
                    if rad <= 0.0 {
                        continue;
                    }
                    rad = rad.sqrt();
                    points.clear();
                    geom_sample_circle(cen, rad, south_dir, density, &mut points).ok()?;
                    if points.is_empty() {
                        continue;
                    }
                    for &point in &points {
                        let mut bail = false;
                        for v in &vectors {
                            let dt2 = point.dot(*v);
                            if dt2 >= 0.0 {
                                bail = true;
                                break;
                            }
                        }
                        if bail {
                            continue;
                        }
                        let point = point + pijk;
                        if hijk < rp && !nears.is_empty() {
                            let mut coll = false;
                            for &np in &nears {
                                let p = &probes[np];
                                if point.distance_squared(p.point) < rp2 {
                                    coll = true;
                                    break;
                                }
                            }
                            if coll {
                                continue;
                            }
                        }
                        let mut mc = 0usize;
                        let mut dmin = 2.0 * rp;
                        for kk in 0..3 {
                            let d = point.distance(atoms[aidx[kk]].coor) - atoms[aidx[kk]].radius;
                            if d < dmin {
                                dmin = d;
                                mc = kk;
                            }
                        }
                        let atom_index = aidx[mc];
                        let molecule = atoms[atom_index].molecule;
                        let pcen = pijk;
                        let outnml = if rp <= 0.0 {
                            point - atoms[atom_index].coor
                        } else {
                            (pcen - point) / rp
                        };
                        let other_mol = u8::from(molecule == 0);
                        let mut buried = false;
                        for b in atoms {
                            if b.molecule != other_mol {
                                continue;
                            }
                            let erl = b.radius + rp;
                            let d = pcen.distance_squared(b.coor);
                            if d <= erl * erl {
                                buried = true;
                                break;
                            }
                        }
                        let dot = Dot {
                            coor: point,
                            outnml,
                            buried,
                        };
                        if molecule == 0 {
                            d0.push(dot);
                        } else {
                            d1.push(dot);
                        }
                    }
                }
                if d0.is_empty() && d1.is_empty() {
                    None
                } else {
                    Some((d0, d1))
                }
            })
            .collect();
        for (mut d0, mut d1) in results {
            self.run.dots[0].append(&mut d0);
            self.run.dots[1].append(&mut d1);
        }
        Ok(())
    }

    fn add_dot(&mut self, molecule: u8, coor: Vec3, pcen: Vec3) {
        let outnml = (pcen - coor) / PROBE_RADIUS;
        let buried = self.run.atoms.iter().any(|b| {
            if b.molecule != molecule {
                let erl = b.radius + PROBE_RADIUS;
                let d = pcen.distance_squared(b.coor);
                d <= erl * erl
            } else {
                false
            }
        });
        let dot = Dot {
            coor,
            outnml,
            buried,
        };
        self.run.dots[usize::from(molecule)].push(dot);
    }
}

// Pure geometry helpers for use in parallel closures
fn geom_sample_arc_segment(
    cen: Vec3,
    rad: f64,
    x: Vec3,
    y: Vec3,
    angle: f64,
    density: f64,
    points: &mut Vec<Vec3>,
) -> Result<(), SurfaceCalculatorError> {
    if rad <= 0.0 {
        points.clear();
        return Ok(());
    }
    let delta = 1.0 / (density.sqrt() * rad);
    if !delta.is_finite() || delta <= 0.0 || !angle.is_finite() || angle < 0.0 {
        return Err(SurfaceCalculatorError::InvalidInput(
            "surface sampling geometry must be finite and non-negative".into(),
        ));
    }
    points.clear();
    let sample_count = (angle / delta + 0.5).floor() as usize;
    points.reserve(sample_count);
    for step in 0..sample_count {
        let a = delta * (step as f64 + 0.5);
        let c = rad * a.cos();
        let s = rad * a.sin();
        points.push(cen + x * c + y * s);
    }
    Ok(())
}

fn geom_sample_arc(
    cen: Vec3,
    rad: f64,
    axis: Vec3,
    density: f64,
    x: Vec3,
    v: Vec3,
    points: &mut Vec<Vec3>,
) -> Result<(), SurfaceCalculatorError> {
    let y = axis.cross(x);
    let dt1 = v.dot(x);
    let dt2 = v.dot(y);
    let mut angle = dt2.atan2(dt1);
    if angle < 0.0 {
        angle += 2.0 * PI;
    }
    geom_sample_arc_segment(cen, rad, x, y, angle, density, points)
}

fn geom_sample_circle(
    cen: Vec3,
    rad: f64,
    axis: Vec3,
    density: f64,
    points: &mut Vec<Vec3>,
) -> Result<(), SurfaceCalculatorError> {
    let mut v1 = Vec3::new(
        axis.y * axis.y + axis.z * axis.z,
        axis.x * axis.x + axis.z * axis.z,
        axis.x * axis.x + axis.y * axis.y,
    );
    v1.normalize();
    let dt = v1.dot(axis);
    if dt.abs() > 0.99 {
        v1 = Vec3::new(1.0, 0.0, 0.0);
    }
    let mut v2 = axis.cross(v1);
    v2.normalize();
    let mut x = axis.cross(v2);
    x.normalize();
    let y = axis.cross(x);
    geom_sample_arc_segment(cen, rad, x, y, 2.0 * PI, density, points)
}
