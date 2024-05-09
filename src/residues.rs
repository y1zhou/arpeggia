use nalgebra as na;
use pdbtbx::*;
use rayon::prelude::*;

#[derive(Clone, Copy, Debug)]
pub struct Ring {
    pub center: na::Vector3<f64>,
    pub normal: na::Vector3<f64>,
}

pub trait ResidueExt {
    /// The residue one-letter code, or `None` if it's not an amino acid.
    fn resn(&self) -> Option<&str>;
    /// Return the atoms in the aromatic ring of the residue.
    fn ring_atoms(&self) -> Vec<&Atom>;
    /// Return the center and normal of the aromatic ring of the residue.
    fn ring_center_and_normal(&self) -> Option<Ring>;
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
        let res_name = self
            .name()
            .unwrap_or("Some conformers have different names!");

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

    fn ring_center_and_normal(&self) -> Option<Ring> {
        let ring_atoms = self.ring_atoms();
        if ring_atoms.is_empty() {
            return None;
        }

        // Construct 3*N matrix of ring atom coordinates
        let mut atom_coords: na::Matrix3xX<f64> = na::Matrix3xX::<f64>::from_iterator(
            ring_atoms.len(),
            ring_atoms.iter().flat_map(|atom| {
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
        // Sanity check the normal is orthogonal to the plane vectors
        // let dot_products = atom_coords.transpose() * normal;
        // debug!("Dot products:\n{:?}", dot_products);

        Some(Ring { center, normal })
    }
}
