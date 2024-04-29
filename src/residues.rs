use pdbtbx::*;
use std::collections::HashSet;

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

pub fn is_hydrogen_acceptor(residue: &Residue, atom: &Atom) -> bool {
    // all the carbonyl oxygens in the main chain and on the terminals
    let oxygens = HashSet::from(["O", "OXT"]);
    if oxygens.contains(atom.name()) && residue.name().unwrap() != "HOH" {
        return true;
    }
    matches!(
        (residue.name().unwrap(), atom.name()),
        ("ASN", "OD1")
            // | ("ASN", "ND2")
            | ("ASP", "OD1")
            | ("ASP", "OD2")
            | ("GLN", "OE1")
            | ("GLU", "OE1")
            | ("GLU", "OE2")
            | ("HIS", "ND1")
            | ("HIS", "NE2")
            | ("SER", "OG")
            | ("THR", "OG1")
            | ("TYR", "OH")
            | ("MET", "SD") // 10.1021/jz300207k and 10.1002/prot.22327
            | ("CYS", "SG") // 10.1002/prot.22327
    )
}

pub fn is_hydrogen_donor(residue: &Residue, atom: &Atom) -> bool {
    // All amide niteogens in the main chain except proline
    if (residue.name().unwrap() != "PRO") && atom.name() == "N" {
        return true;
    }
    matches!(
        (residue.name().unwrap(), atom.name()),
        ("ARG", "NE")
            | ("ARG", "NH1")
            | ("ARG", "NH2")
            | ("ASN", "ND2")
            | ("GLN", "NE2")
            | ("HIS", "ND1")
            | ("HIS", "NE2")
            | ("LYS", "NZ")
            | ("SER", "OG")
            | ("THR", "OG1")
            | ("TRP", "NE1")
            | ("TYR", "OH")
            | ("CYS", "SG") // 10.1002/prot.22327
    )
}
