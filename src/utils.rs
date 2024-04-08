use pdbtbx::*;

pub trait ResidueExt {
    fn resn(&self) -> Option<&str>;
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
            _ => "X",
        };

        match aa_code {
            "X" => None,
            _ => Some(aa_code),
        }
    }
}

pub fn load_model(input_file: &String) -> (PDB, Vec<PDBError>) {
    // Load file as complex structure
    let (mut pdb, errors) = pdbtbx::open(input_file, StrictnessLevel::Loose).unwrap();

    // Remove non-protein residues from model
    pdb.remove_residues_by(|res| res.resn().is_none());

    (pdb, errors)
}
