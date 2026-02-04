//! Connolly molecular surface generator for SC calculations.
//! Ported from <https://github.com/cytokineking/sc-rs>

use std::cmp::Ordering;
use std::f64::consts::PI;
use std::fmt;

use super::atomic_radii::{embedded_atomic_radii, wildcard_match};
use super::settings::Settings;
use super::types::*;
use super::vector3::Vec3;
use rayon::prelude::*;

/// Error type for surface calculation operations.
#[derive(Debug)]
pub enum SurfaceCalculatorError {
    /// No atoms defined
    NoAtoms,
    /// Failed to read radii
    Io(std::io::Error),
    /// Overlapping atoms detected
    Coincident(String),
    /// Sampling limit exceeded
    TooManySubdivisions,
}

impl fmt::Display for SurfaceCalculatorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SurfaceCalculatorError::NoAtoms => write!(f, "No atoms defined"),
            SurfaceCalculatorError::Io(e) => write!(f, "Failed to read radii: {e}"),
            SurfaceCalculatorError::Coincident(msg) => {
                write!(f, "Overlapping atoms detected: {msg}")
            }
            SurfaceCalculatorError::TooManySubdivisions => write!(f, "Sampling limit exceeded"),
        }
    }
}

impl std::error::Error for SurfaceCalculatorError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            SurfaceCalculatorError::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for SurfaceCalculatorError {
    fn from(err: std::io::Error) -> Self {
        SurfaceCalculatorError::Io(err)
    }
}

/// Type alias for parallel neighbor computation results
type NeighborResult = Result<Vec<(Vec<usize>, Vec<usize>, bool)>, SurfaceCalculatorError>;

pub struct SurfaceGenerator {
    pub settings: Settings,
    radii: Vec<AtomRadius>,
    pub(crate) run: RunState,
}

#[derive(Clone, Default)]
pub(crate) struct RunState {
    pub atoms: Vec<ScAtom>,
    pub probes: Vec<Probe>,
    pub dots: [Vec<Dot>; 2],
    pub trimmed_dots: [Vec<usize>; 2],
    pub results: Results,
    pub radmax: f64,
}

impl Default for SurfaceGenerator {
    fn default() -> Self {
        Self::new()
    }
}

impl SurfaceGenerator {
    pub fn new() -> Self {
        Self {
            settings: Settings::default(),
            radii: Vec::new(),
            run: RunState::default(),
        }
    }

    pub fn init(&mut self) {
        if self.radii.is_empty() {
            self.radii = embedded_atomic_radii();
        }
    }

    pub fn get_atom_radius(&self, resn: &str, atomn: &str) -> Result<f64, SurfaceCalculatorError> {
        // First try the embedded radii table
        for radius in &self.radii {
            if !wildcard_match(resn, &radius.residue) {
                continue;
            }
            if !wildcard_match(atomn, &radius.atom) {
                continue;
            }
            return Ok(radius.radius);
        }
        Err(SurfaceCalculatorError::Io(std::io::Error::other(format!(
            "No radius for {resn}:{atomn}"
        ))))
    }

    pub fn assign_attention_numbers(&mut self) {
        self.run.results.surfaces[0].n_buried_atoms = 0;
        self.run.results.surfaces[0].n_blocked_atoms = 0;
        self.run.results.surfaces[1].n_buried_atoms = 0;
        self.run.results.surfaces[1].n_blocked_atoms = 0;

        let sep2 = self.settings.separation_cutoff * self.settings.separation_cutoff;
        let atoms_len = self.run.atoms.len();
        let snapshot: Vec<(usize, usize, f64)> = (0..atoms_len)
            .into_par_iter()
            .map(|i| {
                let a1 = &self.run.atoms[i];
                let mut dist_min2 = f64::INFINITY;
                for j in 0..atoms_len {
                    let a2 = &self.run.atoms[j];
                    if a2.molecule == a1.molecule {
                        continue;
                    }
                    let r2 = a1.distance_squared(a2);
                    if r2 < dist_min2 {
                        dist_min2 = r2;
                    }
                }
                (i, a1.molecule, dist_min2)
            })
            .collect();
        for (i, mol, dist_min2) in snapshot {
            let a1 = &mut self.run.atoms[i];
            if dist_min2 >= sep2 {
                a1.attention = Attention::Far;
                self.run.results.surfaces[mol].n_blocked_atoms += 1;
            } else {
                a1.attention = Attention::Buried;
                self.run.results.surfaces[mol].n_buried_atoms += 1;
            }
        }
    }

    pub(crate) fn generate_molecular_surfaces(&mut self) -> Result<(), SurfaceCalculatorError> {
        if self.run.atoms.is_empty() {
            return Err(SurfaceCalculatorError::NoAtoms);
        }
        self.calc_dots_for_all_atoms()?;
        Ok(())
    }

    fn calc_dots_for_all_atoms(&mut self) -> Result<(), SurfaceCalculatorError> {
        self.run.radmax = 0.0;
        for a in &self.run.atoms {
            if a.radius > self.run.radmax {
                self.run.radmax = a.radius;
            }
        }
        self.compute_neighbors_all_parallel()?;
        let len = self.run.atoms.len();
        for i in 0..len {
            let att = self.run.atoms[i].attention;
            if matches!(att, Attention::Far) {
                continue;
            }
            if matches!(self.run.atoms[i].attention, Attention::Consider)
                && self.run.atoms[i].buried_by_indices.is_empty()
            {
                continue;
            }
            self.build_probes(i)?;
        }
        self.generate_contact_surface_parallel()?;
        if self.settings.rp > 0.0 {
            self.generate_concave_surface_parallel()?;
        }
        Ok(())
    }

    fn compute_neighbors_all_parallel(&mut self) -> Result<(), SurfaceCalculatorError> {
        let rp = self.settings.rp;
        let atoms: &Vec<ScAtom> = &self.run.atoms;
        let atomi_table = atoms
            .iter()
            .enumerate()
            .map(|(i, a)| (a.atomi, i))
            .collect::<std::collections::HashMap<usize, usize>>();
        let results: NeighborResult = atoms
            .par_iter()
            .enumerate()
            .map(|(i, atom1)| {
                let mut neighbor_indices: Vec<usize> = Vec::new();
                let mut buried_by_indices: Vec<usize> = Vec::new();
                // for (j, atom2) in atoms.iter().enumerate() {
                if i == 1 {
                    println!(
                        "Atom {i} atomi{} checking neighbors {:?}",
                        atom1.atomi, atom1.all_neighbors_atomi
                    );
                }
                for (j, atom2_idx) in atom1.all_neighbors_atomi.iter().enumerate() {
                    // if j == i {
                    //     continue;
                    // }
                    // let atom2_idx = &atom2.atomi;
                    let atom2 = &atoms[*atomi_table.get(atom2_idx).unwrap()];
                    if &atom1.atomi == atom2_idx {
                        continue;
                    }
                    let d2 = atom1.distance_squared(atom2);
                    if atom1.molecule == atom2.molecule {
                        if d2 <= 0.0001 {
                            return Err(SurfaceCalculatorError::Coincident(format!(
                                "{}:{}:{} == {}:{}:{}",
                                atom1.atomi,
                                atom1.resn,
                                atom1.atomn,
                                atom2.atomi,
                                atom2.resn,
                                atom2.atomn
                            )));
                        }
                        let bridge = atom1.radius + atom2.radius + 2.0 * rp;
                        if i == 1 {
                            println!(
                                "atom2 atomi{}: d2={:.3}, bridge2={:.3}",
                                atom2.atomi,
                                d2,
                                bridge * bridge
                            );
                        }
                        if d2 < bridge * bridge {
                            neighbor_indices.push(j);
                        }
                    } else {
                        let bridge = atom1.radius + atom2.radius + 2.0 * rp;
                        if d2 < bridge * bridge {
                            buried_by_indices.push(j);
                        }
                    }
                }
                let center = atom1.coor;
                neighbor_indices.sort_unstable_by(|&a1, &a2| {
                    let d1 = atoms[a1].coor.distance_squared(center);
                    let d2 = atoms[a2].coor.distance_squared(center);
                    if d1 < d2 {
                        Ordering::Less
                    } else if d1 > d2 {
                        Ordering::Greater
                    } else {
                        Ordering::Equal
                    }
                });
                let accessible = neighbor_indices.is_empty();
                Ok((neighbor_indices, buried_by_indices, accessible))
            })
            .collect();
        let outs = results?;
        for (i, (neighbors, buried_by, accessible)) in outs.into_iter().enumerate() {
            let a1 = &mut self.run.atoms[i];
            a1.neighbor_indices = neighbors;
            a1.buried_by_indices = buried_by;
            if accessible {
                a1.accessible = true;
            }
        }
        Ok(())
    }

    fn generate_contact_surface_parallel(&mut self) -> Result<(), SurfaceCalculatorError> {
        let rp = self.settings.rp;
        let atoms: &Vec<ScAtom> = &self.run.atoms;
        let results: Vec<(usize, Vec<Dot>, usize)> = (0..atoms.len())
            .into_par_iter()
            .filter_map(|i| {
                let a_i = &atoms[i];
                let att = a_i.attention;
                if matches!(att, Attention::Far) {
                    return None;
                }
                if matches!(att, Attention::Consider) && a_i.buried_by_indices.is_empty() {
                    return None;
                }
                if !a_i.accessible {
                    return None;
                }
                let neighbors = a_i.neighbor_indices.clone();
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
                let cs = geom_sample_arc(
                    o,
                    radius_i,
                    equatorial_vector,
                    a_i.density,
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
                    let ps =
                        geom_sample_circle(cen, rad, north_dir, a_i.density, &mut points).ok()?;
                    if points.is_empty() {
                        continue;
                    }
                    let area = ps * cs;
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
                        let other_mol = usize::from(a_i.molecule == 0);
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
                            area,
                            buried,
                            kind: DotKind::Contact,
                            atom_index: i,
                        });
                    }
                }
                if dots.is_empty() {
                    None
                } else {
                    let n = dots.len();
                    Some((a_i.molecule, dots, n))
                }
            })
            .collect();
        for (mol, mut dots, n) in results {
            self.run.results.dots.convex += n;
            self.run.dots[mol].append(&mut dots);
        }
        Ok(())
    }

    fn build_probes(&mut self, atom_index: usize) -> Result<(), SurfaceCalculatorError> {
        let expanded_radius_i;
        let neighbor_indices: Vec<usize>;
        {
            let atom1 = &self.run.atoms[atom_index];
            expanded_radius_i = atom1.radius + self.settings.rp;
            neighbor_indices = atom1.neighbor_indices.clone();
        }
        for &j in &neighbor_indices {
            let atom2 = &self.run.atoms[j];
            if atom2.atomi <= self.run.atoms[atom_index].atomi {
                continue;
            }
            let expanded_radius_j = atom2.radius + self.settings.rp;
            let dist_ij = self.run.atoms[atom_index].coor.distance(atom2.coor);
            let unit_axis = (atom2.coor - self.run.atoms[atom_index].coor) / dist_ij;
            let asymmetry_term = (expanded_radius_i * expanded_radius_i
                - expanded_radius_j * expanded_radius_j)
                / dist_ij;
            let midplane_center = (self.run.atoms[atom_index].coor + atom2.coor) * 0.5
                + (unit_axis * (asymmetry_term * 0.5));
            let mut far_term = (expanded_radius_i + expanded_radius_j)
                * (expanded_radius_i + expanded_radius_j)
                - dist_ij * dist_ij;
            if far_term <= 0.0 {
                continue;
            }
            far_term = far_term.sqrt();
            let mut contain_term =
                dist_ij * dist_ij - (self.run.atoms[atom_index].radius - atom2.radius).powi(2);
            if contain_term <= 0.0 {
                continue;
            }
            contain_term = contain_term.sqrt();
            let ring_radius = 0.5 * far_term * contain_term / dist_ij;
            if neighbor_indices.len() <= 1 {
                self.run.atoms[atom_index].accessible = true;
                self.run.atoms[j].accessible = true;
                break;
            }
            self.build_probe_triplets(atom_index, j, unit_axis, midplane_center, ring_radius)?;
            let has_point_cusp = asymmetry_term.abs() < dist_ij;
            let atom2_att = self.run.atoms[j].attention;
            if !matches!(self.run.atoms[atom_index].attention, Attention::Far)
                || (!matches!(atom2_att, Attention::Far) && self.settings.rp > 0.0)
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
        let neighbor_indices = self.run.atoms[atom1_index].neighbor_indices.clone();
        let expanded_radius_i = self.run.atoms[atom1_index].radius + self.settings.rp;
        let atom2 = &self.run.atoms[atom2_index];
        let expanded_radius_j = atom2.radius + self.settings.rp;
        let atom2_natom = atom2.atomi;
        let atom2_coor = atom2.coor;
        let atom2_att = atom2.attention;
        let mut made_probe = false;
        for &k in &neighbor_indices {
            let atom3 = &self.run.atoms[k];
            if atom3.atomi <= atom2_natom {
                continue;
            }
            let expanded_radius_k = atom3.radius + self.settings.rp;
            let dist_jk = atom2_coor.distance(atom3.coor);
            if dist_jk >= expanded_radius_j + expanded_radius_k {
                continue;
            }
            let dist_ik = self.run.atoms[atom1_index].coor.distance(atom3.coor);
            if dist_ik >= expanded_radius_i + expanded_radius_k {
                continue;
            }
            if matches!(self.run.atoms[atom1_index].attention, Attention::Far)
                && matches!(atom2_att, Attention::Far)
                && matches!(atom3.attention, Attention::Far)
            {
                continue;
            }
            let unit_axis_ik = (atom3.coor - self.run.atoms[atom1_index].coor) / dist_ik;
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
            let midpoint_ik = (self.run.atoms[atom1_index].coor + atom3.coor) * 0.5
                + unit_axis_ik * (asymmetry_term_ik * 0.5);
            let mut componentwise = midpoint_ik - midplane_center;
            componentwise = Vec3::new(
                unit_axis_ik.x * componentwise.x,
                unit_axis_ik.y * componentwise.y,
                unit_axis_ik.z * componentwise.z,
            );
            let component_sum = componentwise.x + componentwise.y + componentwise.z;
            let torus_center = midplane_center + perp_tangent * (component_sum / sin_wedge);
            let mut height = expanded_radius_i * expanded_radius_i
                - torus_center.distance_squared(self.run.atoms[atom1_index].coor);
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
        let density = f64::midpoint(
            self.run.atoms[atom1_index].density,
            self.run.atoms[atom2_index].density,
        );
        let expanded_radius_i = self.run.atoms[atom1_index].radius + self.settings.rp;
        let expanded_radius_j = self.run.atoms[atom2_index].radius + self.settings.rp;
        let roll_circle_radius_i =
            ring_radius * self.run.atoms[atom1_index].radius / expanded_radius_i;
        let roll_circle_radius_j =
            ring_radius * self.run.atoms[atom2_index].radius / expanded_radius_j;
        let mut belt_radius = ring_radius - self.settings.rp;
        if belt_radius <= 0.0 {
            belt_radius = 0.0;
        }
        let mean_radius = (roll_circle_radius_i + 2.0 * belt_radius + roll_circle_radius_j) / 4.0;
        let eccentricity = mean_radius / ring_radius;
        let effective_density = eccentricity * eccentricity * density;
        let mut subs: Vec<Vec3> = Vec::new();
        let ts = self.sample_circle(
            midplane_center,
            ring_radius,
            unit_axis,
            effective_density,
            &mut subs,
        )?;
        if subs.is_empty() {
            return Ok(());
        }
        let atom2_natom = self.run.atoms[atom2_index].atomi;
        for sub in subs {
            let mut tooclose = false;
            for &ni in &neighbors {
                let neighbor = &self.run.atoms[ni];
                if neighbor.atomi == atom2_natom {
                    continue;
                }
                let expanded_neighbor_radius = neighbor.radius + self.settings.rp;
                let d2 = sub.distance_squared(neighbor.coor);
                if d2 < expanded_neighbor_radius * expanded_neighbor_radius {
                    tooclose = true;
                    break;
                }
            }
            if tooclose {
                continue;
            }
            let ring_point = sub;
            self.run.atoms[atom1_index].accessible = true;
            self.run.atoms[atom2_index].accessible = true;
            let vec_pi = (self.run.atoms[atom1_index].coor - ring_point) / expanded_radius_i;
            let vec_pj = (self.run.atoms[atom2_index].coor - ring_point) / expanded_radius_j;
            let mut toroid_axis = vec_pi.cross(vec_pj);
            toroid_axis.normalize();
            let mut cusp_term = self.settings.rp * self.settings.rp - ring_radius * ring_radius;
            let has_cusp_point = cusp_term > 0.0 && has_point_cusp;
            let (arc_end_i, arc_end_j) = if has_cusp_point {
                cusp_term = cusp_term.sqrt();
                let qij = midplane_center - unit_axis * cusp_term;
                (((qij - ring_point) / self.settings.rp), Vec3::zero())
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
                let ps = self.sample_arc(
                    ring_point,
                    self.settings.rp,
                    toroid_axis,
                    density,
                    vec_pi,
                    arc_end_i,
                    &mut points,
                )?;
                for &point in &points {
                    let area = ps * ts * distance_point_to_line(midplane_center, unit_axis, point)
                        / ring_radius;
                    self.run.results.dots.toroidal += 1;
                    let molecule = self.run.atoms[atom1_index].molecule;
                    self.add_dot(
                        molecule,
                        DotKind::Reentrant,
                        point,
                        area,
                        ring_point,
                        atom1_index,
                    );
                }
            }
            let atom2_attention = self.run.atoms[atom2_index].attention;
            if !matches!(atom2_attention, Attention::Far) {
                let mut points: Vec<Vec3> = Vec::new();
                let ps = self.sample_arc(
                    ring_point,
                    self.settings.rp,
                    toroid_axis,
                    density,
                    arc_end_j,
                    vec_pj,
                    &mut points,
                )?;
                for &point in &points {
                    let area = ps * ts * distance_point_to_line(midplane_center, unit_axis, point)
                        / ring_radius;
                    self.run.results.dots.toroidal += 1;
                    let molecule2 = self.run.atoms[atom2_index].molecule;
                    self.add_dot(
                        molecule2,
                        DotKind::Reentrant,
                        point,
                        area,
                        ring_point,
                        atom2_index,
                    );
                }
            }
        }
        Ok(())
    }

    fn check_atom_collision2_idx(
        &self,
        probe_center: Vec3,
        atom1_index: usize,
        atom2_index: usize,
        neighbor_indices: &Vec<usize>,
    ) -> bool {
        let atom1_natom = self.run.atoms[atom1_index].atomi;
        let atom2_natom = self.run.atoms[atom2_index].atomi;
        for &ni in neighbor_indices {
            let neighbor = &self.run.atoms[ni];
            if neighbor.atomi == atom1_natom || neighbor.atomi == atom2_natom {
                continue;
            }
            if probe_center.distance_squared(neighbor.coor)
                <= (neighbor.radius + self.settings.rp).powi(2)
            {
                return true;
            }
        }
        false
    }

    fn generate_concave_surface_parallel(&mut self) -> Result<(), SurfaceCalculatorError> {
        let rp = self.settings.rp;
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
        let results: Vec<(Vec<Dot>, Vec<Dot>, usize)> = (0..probes.len())
            .into_par_iter()
            .filter_map(|i| {
                let probe = &probes[i];
                let aidx = probe.atom_indices;
                if matches!(atoms[aidx[0]].attention, Attention::Consider)
                    && matches!(atoms[aidx[1]].attention, Attention::Consider)
                    && matches!(atoms[aidx[2]].attention, Attention::Consider)
                {
                    return None;
                }
                let pijk = probe.point;
                let uijk = probe.alt;
                let hijk = probe.height;
                let density =
                    (atoms[aidx[0]].density + atoms[aidx[1]].density + atoms[aidx[2]].density)
                        / 3.0;
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
                let cs =
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
                    let ps = geom_sample_circle(cen, rad, south_dir, density, &mut points).ok()?;
                    if points.is_empty() {
                        continue;
                    }
                    let area = ps * cs;
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
                        let other_mol = usize::from(molecule == 0);
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
                            area,
                            buried,
                            kind: DotKind::Cavity,
                            atom_index,
                        };
                        if molecule == 0 {
                            d0.push(dot);
                        } else {
                            d1.push(dot);
                        }
                    }
                }
                let n = d0.len() + d1.len();
                if n == 0 { None } else { Some((d0, d1, n)) }
            })
            .collect();
        for (mut d0, mut d1, n) in results {
            self.run.results.dots.concave += n;
            self.run.dots[0].append(&mut d0);
            self.run.dots[1].append(&mut d1);
        }
        Ok(())
    }

    fn add_dot(
        &mut self,
        molecule: usize,
        kind: DotKind,
        coor: Vec3,
        area: f64,
        pcen: Vec3,
        atom_index: usize,
    ) {
        let atom = &self.run.atoms[atom_index];
        let outnml = if self.settings.rp <= 0.0 {
            coor - atom.coor
        } else {
            (pcen - coor) / self.settings.rp
        };
        let mut buried = false;
        let other_mol = usize::from(molecule == 0);
        for b in &self.run.atoms {
            if b.molecule != other_mol {
                continue;
            }
            let erl = b.radius + self.settings.rp;
            let d = pcen.distance_squared(b.coor);
            if d <= erl * erl {
                buried = true;
                break;
            }
        }
        let dot = Dot {
            coor,
            outnml,
            area,
            buried,
            kind,
            atom_index,
        };
        self.run.dots[molecule].push(dot);
    }

    #[allow(clippy::too_many_arguments)]
    fn sample_arc(
        &self,
        cen: Vec3,
        rad: f64,
        axis: Vec3,
        density: f64,
        x: Vec3,
        v: Vec3,
        points: &mut Vec<Vec3>,
    ) -> Result<f64, SurfaceCalculatorError> {
        let y = axis.cross(x);
        let dt1 = v.dot(x);
        let dt2 = v.dot(y);
        let mut angle = dt2.atan2(dt1);
        if angle < 0.0 {
            angle += 2.0 * PI;
        }
        sample_arc_segment(cen, rad, x, y, angle, density, points)
    }

    fn sample_circle(
        &self,
        cen: Vec3,
        rad: f64,
        axis: Vec3,
        density: f64,
        points: &mut Vec<Vec3>,
    ) -> Result<f64, SurfaceCalculatorError> {
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
        sample_arc_segment(cen, rad, x, y, 2.0 * PI, density, points)
    }
}

fn distance_point_to_line(cen: Vec3, axis: Vec3, pnt: Vec3) -> f64 {
    let vec = pnt - cen;
    let dt = vec.dot(axis);
    let mut d2 = vec.magnitude_squared() - dt * dt;
    if d2 < 0.0 {
        d2 = 0.0;
    }
    d2.sqrt()
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
) -> Result<f64, SurfaceCalculatorError> {
    if rad <= 0.0 {
        points.clear();
        return Ok(0.0);
    }
    let delta = 1.0 / (density.sqrt() * rad);
    let mut a = -delta / 2.0;
    points.clear();
    for _ in 0..100000 {
        a += delta;
        if a > angle {
            break;
        }
        let c = rad * a.cos();
        let s = rad * a.sin();
        points.push(cen + x * c + y * s);
    }
    if a + delta < angle {
        return Err(SurfaceCalculatorError::TooManySubdivisions);
    }
    let ps = if points.is_empty() {
        0.0
    } else {
        rad * angle / (points.len() as f64)
    };
    Ok(ps)
}

fn geom_sample_arc(
    cen: Vec3,
    rad: f64,
    axis: Vec3,
    density: f64,
    x: Vec3,
    v: Vec3,
    points: &mut Vec<Vec3>,
) -> Result<f64, SurfaceCalculatorError> {
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
) -> Result<f64, SurfaceCalculatorError> {
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

#[allow(clippy::too_many_arguments)]
fn sample_arc_segment(
    cen: Vec3,
    rad: f64,
    x: Vec3,
    y: Vec3,
    angle: f64,
    density: f64,
    points: &mut Vec<Vec3>,
) -> Result<f64, SurfaceCalculatorError> {
    if rad <= 0.0 {
        points.clear();
        return Ok(0.0);
    }
    let delta = 1.0 / (density.sqrt() * rad);
    let mut a = -delta / 2.0;
    points.clear();
    for _ in 0..100000 {
        a += delta;
        if a > angle {
            break;
        }
        let c = rad * a.cos();
        let s = rad * a.sin();
        points.push(cen + x * c + y * s);
    }
    if a + delta < angle {
        return Err(SurfaceCalculatorError::TooManySubdivisions);
    }
    let ps = if points.is_empty() {
        0.0
    } else {
        rad * angle / (points.len() as f64)
    };
    Ok(ps)
}
