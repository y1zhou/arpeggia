use std::collections::HashSet;

use pdbtbx::*;

pub fn is_hydrogen_acceptor(res_name: &str, atom_name: &str) -> bool {
    // all the carbonyl oxygens in the main chain and on the terminals
    let oxygens = HashSet::from(["O", "OXT"]);
    if oxygens.contains(atom_name) && res_name != "HOH" {
        return true;
    }
    matches!(
        (res_name, atom_name),
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

pub fn is_hydrogen_donor(res_name: &str, atom_name: &str) -> bool {
    // All amide niteogens in the main chain except proline
    if atom_name == "CA" {
        return true;
    }
    matches!(
        (res_name, atom_name),
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
