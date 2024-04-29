use std::collections::HashSet;

use pdbtbx::*;

pub fn is_hydrogen_acceptor(rotamer: &Conformer, atom: &Atom) -> bool {
    // all the carbonyl oxygens in the main chain and on the terminals
    let oxygens = HashSet::from(["O", "OXT"]);
    if oxygens.contains(atom.name()) && rotamer.name() != "HOH" {
        return true;
    }
    matches!(
        (rotamer.name(), atom.name()),
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

pub fn is_hydrogen_donor(rotamer: &Conformer, atom: &Atom) -> bool {
    // All amide niteogens in the main chain except proline
    if atom.name() == "CA" {
        return true;
    }
    matches!(
        (rotamer.name(), atom.name()),
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

pub fn is_weak_hydrogen_donor(atom: &Atom) -> bool {
    // All the non-carbonyl carbon atoms
    (atom.element().unwrap().to_string() == "C") && atom.name() != "C"
}
