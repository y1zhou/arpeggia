use nalgebra as na;
use pdbtbx::*;
use rayon::prelude::*;

/// The struct for a residue identifier
#[derive(Debug, Hash, PartialEq, Eq, Clone)]
pub struct ResidueId<'a> {
    /// Model identifier
    pub model: usize,
    /// Chain identifier
    pub chain: &'a str,
    /// Residue index
    pub resi: isize,
    /// Residue insertion code
    pub insertion: &'a str,
    /// Alternate location identifier
    pub altloc: &'a str,
    /// Residue name
    pub resn: &'a str,
}

/// The struct for a plane in 3D space
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Plane {
    pub center: na::Vector3<f64>,
    pub normal: na::Vector3<f64>,
}

impl Plane {
    /// Calculate the distance from a point to the plane center
    pub fn point_vec_dist(&self, point: &na::Vector3<f64>) -> f64 {
        (point - self.center).norm()
    }

    pub fn point_dist(&self, point: &(f64, f64, f64)) -> f64 {
        let atom_point = na::Vector3::new(point.0, point.1, point.2);
        self.point_vec_dist(&atom_point)
    }

    /// Calculate the angle between the plane normal and the vector pointing from
    /// the plane center to the point.
    pub fn point_vec_angle(&self, point: &na::Vector3<f64>) -> f64 {
        let v = point - self.center;
        let mut rad = (self.normal.dot(&v) / (self.normal.norm() * v.norm())).acos();

        // Convert to degrees
        if rad > std::f64::consts::FRAC_PI_2 {
            rad = std::f64::consts::PI - rad;
        }
        rad.to_degrees()
    }

    #[allow(dead_code)]
    pub fn point_angle(&self, point: &(f64, f64, f64)) -> f64 {
        let atom_point = na::Vector3::new(point.0, point.1, point.2);
        self.point_vec_angle(&atom_point)
    }

    /// Calculate the angle between two planes
    pub fn dihedral(&self, plane: &Plane) -> f64 {
        let mut rad =
            (self.normal.dot(&plane.normal) / (self.normal.norm() * plane.normal.norm())).acos();

        // Convert to degrees
        if rad > std::f64::consts::FRAC_PI_2 {
            rad = std::f64::consts::PI - rad;
        }
        rad.to_degrees()
    }
}

impl<'a> ResidueId<'a> {
    pub fn new(
        model: usize,
        chain: &'a str,
        resi: isize,
        insertion: &'a str,
        altloc: &'a str,
        resn: &'a str,
    ) -> Self {
        Self {
            model,
            chain,
            resi,
            insertion,
            altloc,
            resn,
        }
    }

    /// Helper function to convert an [`pdbtbx::AtomConformerResidueChainModel`] to a residue identifier
    pub fn from_hier(hier: &'a AtomConformerResidueChainModel) -> Self {
        let (resi, insertion) = hier.residue().id();
        let altloc = hier.conformer().alternative_location();
        Self::new(
            hier.model().serial_number(),
            hier.chain().id(),
            resi,
            insertion.unwrap_or(""),
            altloc.unwrap_or(""),
            hier.residue().name().unwrap_or(""),
        )
    }
}

pub trait ResidueExt {
    /// The residue one-letter code, or `None` if it's not an amino acid.
    fn resn(&self) -> Option<&str>;

    /// Return the atoms in the aromatic ring of the residue.
    fn ring_atoms(&self) -> Vec<&Atom>;

    /// Return the atoms that form a plane in the side chain.
    /// See page 27 of:
    /// https://cdn.rcsb.org/wwpdb/docs/documentation/file-format/PDB_format_1992.pdf
    fn sc_plane_atoms(&self) -> Vec<&Atom>;

    /// Return the center and normal of the given atoms. If no atoms are given,
    /// all side chain plane atoms are used.
    fn center_and_normal(&self, atoms: Option<Vec<&Atom>>) -> Option<Plane>;
}

impl ResidueExt for Residue {
    fn resn(&self) -> Option<&str> {
        let aa_code = match self.name().unwrap().to_uppercase().as_str() {
            "ALA" => "A",
            "ARG" => "R",
            "ASN" => "N",
            "ASP" => "D",
            "CYS" => "C",
            "GLN" => "Q",
            "GLU" => "E",
            "GLY" => "G",
            "HIS" => "H",
            "ILE" => "I",
            "LEU" => "L",
            "LYS" => "K",
            "MET" => "M",
            "PHE" => "F",
            "PRO" => "P",
            "SER" => "S",
            "THR" => "T",
            "TRP" => "W",
            "TYR" => "Y",
            "VAL" => "V",
            "HOH" => "O", // water
            _ => "X",
        };

        match aa_code {
            "X" => None,
            _ => Some(aa_code),
        }
    }

    fn ring_atoms(&self) -> Vec<&Atom> {
        let res_name = self.name().unwrap_or(""); // Some conformers have different names

        match res_name {
            "HIS" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CG" | "ND1" | "CE1" | "NE2" | "CD2"))
                .collect(),
            "PHE" | "TYR" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CG" | "CD1" | "CD2" | "CE1" | "CE2" | "CZ"))
                .collect(),
            "TRP" => self
                .par_atoms()
                .filter(|atom| {
                    matches!(
                        atom.name(),
                        "CG" | "CD1" | "CD2" | "NE1" | "CE2" | "CE3" | "CZ2" | "CZ3" | "CH2"
                    )
                })
                .collect(),
            _ => vec![],
        }
    }

    fn sc_plane_atoms(&self) -> Vec<&Atom> {
        let res_name = self.name().unwrap_or("");

        match res_name {
            "ARG" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "NE" | "CZ" | "NH1" | "NH2"))
                .collect(),
            "ASN" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CB" | "CG" | "OD1" | "ND2"))
                .collect(),
            "ASP" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CB" | "CG" | "OD1" | "OD2"))
                .collect(),
            "CYS" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CA" | "CB" | "SG"))
                .collect(),
            "GLU" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CG" | "CD" | "OE1" | "OE2"))
                .collect(),
            "GLN" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CG" | "CD" | "OE1" | "NE2"))
                .collect(),
            "HIS" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CG" | "ND1" | "CE1" | "NE2" | "CD2"))
                .collect(),
            "ILE" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CB" | "CG1" | "CG2" | "CD1"))
                .collect(),
            "LEU" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CB" | "CG" | "CD1" | "CD2"))
                .collect(),
            "LYS" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CG" | "CD" | "CE" | "NZ"))
                .collect(),
            "MET" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CG" | "SD" | "CE"))
                .collect(),
            "PHE" | "TYR" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CG" | "CD1" | "CD2" | "CE1" | "CE2" | "CZ"))
                .collect(),
            "PRO" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "N" | "CA" | "CB" | "CG" | "CD"))
                .collect(),
            "SER" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CA" | "CB" | "OG"))
                .collect(),
            "THR" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CA" | "CB" | "OG1" | "CG2"))
                .collect(),
            "TRP" => self
                .par_atoms()
                .filter(|atom| {
                    matches!(
                        atom.name(),
                        "CG" | "CD1" | "CD2" | "NE1" | "CE2" | "CE3" | "CZ2" | "CZ3" | "CH2"
                    )
                })
                .collect(),
            "VAL" => self
                .par_atoms()
                .filter(|atom| matches!(atom.name(), "CA" | "CB" | "CG1" | "CG2"))
                .collect(),
            // Nothing for alanine and glycine
            _ => vec![],
        }
    }

    fn center_and_normal(&self, atoms: Option<Vec<&Atom>>) -> Option<Plane> {
        let sc_atoms = atoms.unwrap_or(self.sc_plane_atoms());

        if sc_atoms.len() < 3 {
            return None;
        }

        // Construct 3*N matrix of ring atom coordinates
        let mut atom_coords: na::Matrix3xX<f64> = na::Matrix3xX::<f64>::from_iterator(
            sc_atoms.len(),
            sc_atoms.iter().flat_map(|atom| {
                let coord = atom.pos();
                [coord.0, coord.1, coord.2].into_iter()
            }),
        );
        let center = atom_coords.column_mean();

        // Center the matrix and perform SVD
        // Ref: https://en.wikipedia.org/wiki/Singular_value_decomposition#Total_least_squares_minimization
        // Ref: https://math.stackexchange.com/a/2810220
        for i in 0..atom_coords.ncols() {
            atom_coords.set_column(i, &(atom_coords.column(i) - center));
        }

        let svd = atom_coords.svd(true, true);
        let normal = svd.u.unwrap().column(2).clone_owned();

        Some(Plane { center, normal })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::load_model;

    #[test]
    fn test_plane_dist_angles() {
        let plane_x = Plane {
            center: na::Vector3::new(0.0, 0.0, 0.0),
            normal: na::Vector3::new(0.0, 0.0, 1.0),
        };
        let point = na::Vector3::new(0.0, 1.0, 1.0);

        // Test a plane parallel to the x-y plane
        let parallel_x = Plane {
            center: point,
            normal: na::Vector3::new(0.0, 0.0, -1.0),
        };
        assert!((plane_x.point_vec_dist(&point) - 2.0_f64.sqrt()).abs() < 1e-6);
        assert!((plane_x.point_vec_angle(&point) - 45.0) < 1e-6);
        assert!((parallel_x.point_vec_angle(&plane_x.center) - 45.0) < 1e-6);
        assert!(plane_x.dihedral(&parallel_x) < 1e-6);

        // Test a plane perpendicular to the x-y plane
        let perpendicular_x = Plane {
            center: point,
            normal: na::Vector3::new(1.0, 0.0, 0.0),
        };
        assert!((plane_x.point_vec_angle(&point) - 90.0) < 1e-6);
        assert!((perpendicular_x.point_vec_angle(&plane_x.center) - 90.0) < 1e-6);
        assert!((plane_x.dihedral(&perpendicular_x) - 90.0) < 1e-6);
    }

    #[test]
    fn test_residue_ext() {
        let root = env!("CARGO_MANIFEST_DIR");
        let path = format!("{}/{}", root, "test-data/1ubq.pdb");

        let (pdb, _) = load_model(&path);

        // First Met has no rings
        let residue = pdb.residues().next().unwrap();
        assert_eq!(residue.resn(), Some("M"));
        assert_eq!(residue.ring_atoms().len(), 0);
        assert_eq!(residue.sc_plane_atoms().len(), 3);

        // Phe 4 has a ring
        let res_f4 = pdb.residues().find(|res| res.resn() == Some("F")).unwrap();
        assert_eq!(res_f4.serial_number(), 4);
        let ring_atoms = res_f4.ring_atoms();
        assert_eq!(ring_atoms.len(), 6);
        let f4_ring = res_f4.center_and_normal(Some(ring_atoms.clone())).unwrap();
        let f4_ring_default = res_f4.center_and_normal(None).unwrap();
        assert!(f4_ring == f4_ring_default);
        assert!(
            f4_ring.center.relative_eq(
                &na::Vector3::new(24.96883333, 34.687, 6.16233333),
                f64::EPSILON,
                1e-6
            ),
            "Ring center: {:?}",
            f4_ring.center
        );
        assert!(
            f4_ring.normal.relative_eq(
                &na::Vector3::new(0.53253994, -0.82736044, -0.17853828),
                f64::EPSILON,
                1e-6
            ),
            "Normal vector: {:?}",
            f4_ring.normal
        );
        // Sanity check the normal is orthogonal to the plane vectors
        let atom_coords: na::Matrix3xX<f64> = na::Matrix3xX::<f64>::from_iterator(
            ring_atoms.len(),
            ring_atoms.iter().flat_map(|atom| {
                let coord = atom.pos();
                [
                    coord.0 - f4_ring.center[0],
                    coord.1 - f4_ring.center[1],
                    coord.2 - f4_ring.center[2],
                ]
                .into_iter()
            }),
        );
        let dot_products = atom_coords.transpose() * f4_ring.normal;
        assert!(
            dot_products.abs().mean().abs() < 0.02,
            "Dot products: {:?}\nCenter: {:?}\nNormal: {:?}\nAtoms: {:?}",
            dot_products,
            f4_ring.center,
            f4_ring.normal,
            atom_coords
        );
    }
}
