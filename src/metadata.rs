//! Narrow bond and declared-sequence parsing missing from `pdbtbx`.

use crate::{Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::Path;

type CifLoop = (Vec<String>, Vec<Vec<String>>);

/// An atom identifier used by explicit input connectivity.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct BondEndpoint {
    /// Chain identifier.
    pub chain: String,
    /// Residue sequence number.
    pub residue: isize,
    /// Optional insertion code.
    pub insertion: String,
    /// Atom name.
    pub atom: String,
}

impl BondEndpoint {
    /// Construct an endpoint from a parsed atom identity.
    pub fn from_parts(chain: &str, residue: isize, insertion: &str, atom: &str) -> Self {
        Self {
            chain: clean(chain).to_string(),
            residue,
            insertion: clean_unknown(insertion).to_string(),
            atom: clean(atom).to_string(),
        }
    }

    fn new(chain: &str, residue: &str, insertion: &str, atom: &str) -> Option<Self> {
        Some(Self::from_parts(
            chain,
            clean(residue).parse().ok()?,
            insertion,
            atom,
        ))
    }
}

/// Explicit connectivity and declared polymer sequences from one input file.
#[derive(Clone, Debug, Default)]
pub struct StructureMetadata {
    bonds: HashSet<(BondEndpoint, BondEndpoint)>,
    declared_sequences: Vec<(String, String)>,
}

impl StructureMetadata {
    /// Whether the input explicitly connects two atoms.
    pub fn has_bond(&self, first: &BondEndpoint, second: &BondEndpoint) -> bool {
        self.bonds
            .contains(&ordered_bond(first.clone(), second.clone()))
    }

    /// Chain-oriented declared polymer sequences.
    pub fn declared_sequences(&self) -> &[(String, String)] {
        &self.declared_sequences
    }
}

/// Read only connectivity and polymer-sequence metadata from PDB or mmCIF.
pub fn read_metadata(path: impl AsRef<Path>) -> ArpeggiaResult<Analysis<StructureMetadata>> {
    let path = path.as_ref();
    let input = std::fs::read_to_string(path)?;
    let extension = path
        .extension()
        .and_then(|value| value.to_str())
        .map(str::to_ascii_lowercase);
    let is_cif = matches!(extension.as_deref(), Some("cif" | "mmcif"))
        || input.trim_start().starts_with("data_");
    let (metadata, warnings) = if is_cif {
        parse_mmcif(&input)?
    } else {
        parse_pdb(&input)
    };
    Ok(Analysis::new(metadata, warnings))
}

/// Return declared PDB `SEQRES` or mmCIF entity-polymer sequences by chain.
pub fn get_seqres(path: impl AsRef<Path>) -> ArpeggiaResult<Analysis<Vec<(String, String)>>> {
    Ok(read_metadata(path)?.map(|metadata| metadata.declared_sequences))
}

fn parse_pdb(input: &str) -> (StructureMetadata, Vec<AnalysisWarning>) {
    let mut bonds = HashSet::new();
    let mut atoms = HashMap::new();
    let mut conect = Vec::new();
    let mut sequences = Vec::<(String, Vec<String>)>::new();

    for line in input.lines() {
        match field(line, 0, 6).trim() {
            "ATOM" | "HETATM" => {
                if let (Ok(serial), Some(endpoint)) = (
                    field(line, 6, 11).trim().parse::<usize>(),
                    BondEndpoint::new(
                        field(line, 21, 22),
                        field(line, 22, 26),
                        field(line, 26, 27),
                        field(line, 12, 16),
                    ),
                ) {
                    atoms.insert(serial, endpoint);
                }
            }
            "SSBOND" => {
                if let (Some(first), Some(second)) = (
                    BondEndpoint::new(
                        field(line, 15, 16),
                        field(line, 17, 21),
                        field(line, 21, 22),
                        "SG",
                    ),
                    BondEndpoint::new(
                        field(line, 29, 30),
                        field(line, 31, 35),
                        field(line, 35, 36),
                        "SG",
                    ),
                ) {
                    bonds.insert(ordered_bond(first, second));
                }
            }
            "LINK" => {
                if let (Some(first), Some(second)) = (
                    BondEndpoint::new(
                        field(line, 21, 22),
                        field(line, 22, 26),
                        field(line, 26, 27),
                        field(line, 12, 16),
                    ),
                    BondEndpoint::new(
                        field(line, 51, 52),
                        field(line, 52, 56),
                        field(line, 56, 57),
                        field(line, 42, 46),
                    ),
                ) {
                    bonds.insert(ordered_bond(first, second));
                }
            }
            "CONECT" => {
                let serials = field(line, 6, line.len())
                    .as_bytes()
                    .chunks(5)
                    .filter_map(|chunk| std::str::from_utf8(chunk).ok()?.trim().parse().ok())
                    .collect::<Vec<usize>>();
                if let Some((&first, rest)) = serials.split_first() {
                    conect.extend(rest.iter().map(|&second| (first, second)));
                }
            }
            "SEQRES" => {
                let chain = field(line, 11, 12).trim().to_string();
                let monomers = field(line, 19, line.len())
                    .split_whitespace()
                    .map(str::to_string);
                if let Some((_, existing)) = sequences.iter_mut().find(|(id, _)| id == &chain) {
                    existing.extend(monomers);
                } else {
                    sequences.push((chain, monomers.collect()));
                }
            }
            _ => {}
        }
    }
    for (first, second) in conect {
        if let (Some(first), Some(second)) = (atoms.get(&first), atoms.get(&second)) {
            bonds.insert(ordered_bond(first.clone(), second.clone()));
        }
    }
    metadata_with_sequences(bonds, sequences)
}

fn parse_mmcif(input: &str) -> ArpeggiaResult<(StructureMetadata, Vec<AnalysisWarning>)> {
    let loops = cif_loops(input)?;
    let mut bonds = HashSet::new();
    let mut entities = BTreeMap::<String, BTreeMap<usize, String>>::new();
    let mut chain_entities = Vec::<(String, String)>::new();

    for (tags, rows) in loops {
        let columns: HashMap<&str, usize> = tags
            .iter()
            .enumerate()
            .map(|(index, tag)| (tag.as_str(), index))
            .collect();
        if tags.iter().any(|tag| tag.starts_with("_struct_conn.")) {
            for row in &rows {
                let kind = value(row, &columns, "_struct_conn.conn_type_id").unwrap_or("");
                if !(kind.starts_with("covale") || kind.starts_with("disulf")) {
                    continue;
                }
                if let (Some(first), Some(second)) = (
                    cif_endpoint(row, &columns, 1),
                    cif_endpoint(row, &columns, 2),
                ) {
                    bonds.insert(ordered_bond(first, second));
                }
            }
        } else if tags.iter().any(|tag| tag.starts_with("_entity_poly_seq.")) {
            for row in &rows {
                if let (Some(entity), Some(number), Some(monomer)) = (
                    value(row, &columns, "_entity_poly_seq.entity_id"),
                    value(row, &columns, "_entity_poly_seq.num")
                        .and_then(|number| number.parse::<usize>().ok()),
                    value(row, &columns, "_entity_poly_seq.mon_id"),
                ) {
                    entities
                        .entry(entity.to_string())
                        .or_default()
                        .insert(number, monomer.to_string());
                }
            }
        } else if tags
            .iter()
            .any(|tag| tag.starts_with("_pdbx_poly_seq_scheme."))
        {
            for row in &rows {
                if let (Some(chain), Some(entity)) = (
                    value(row, &columns, "_pdbx_poly_seq_scheme.pdb_strand_id")
                        .or_else(|| value(row, &columns, "_pdbx_poly_seq_scheme.asym_id")),
                    value(row, &columns, "_pdbx_poly_seq_scheme.entity_id"),
                ) {
                    for chain in chain.split(',') {
                        let chain = clean(chain).to_string();
                        if !chain_entities.iter().any(|(known, _)| known == &chain) {
                            chain_entities.push((chain, entity.to_string()));
                        }
                    }
                }
            }
        }
    }

    let sequences = chain_entities
        .into_iter()
        .filter_map(|(chain, entity)| {
            let monomers = entities.get(&entity)?;
            Some((chain, monomers.values().cloned().collect::<Vec<_>>()))
        })
        .collect();
    Ok(metadata_with_sequences(bonds, sequences))
}

fn metadata_with_sequences(
    bonds: HashSet<(BondEndpoint, BondEndpoint)>,
    sequences: Vec<(String, Vec<String>)>,
) -> (StructureMetadata, Vec<AnalysisWarning>) {
    let mut unsupported = BTreeSet::new();
    let declared_sequences = sequences
        .into_iter()
        .map(|(chain, monomers)| {
            let sequence = monomers
                .iter()
                .map(|name| match one_letter(name) {
                    Some(letter) => letter,
                    None => {
                        unsupported.insert(name.trim().to_ascii_uppercase());
                        'X'
                    }
                })
                .collect();
            (chain, sequence)
        })
        .collect();
    let warnings = if unsupported.is_empty() {
        Vec::new()
    } else {
        vec![AnalysisWarning::new(
            WarningCode::UnsupportedMonomer,
            format!(
                "unrecognized declared polymer monomers were mapped to X: {}",
                unsupported.into_iter().collect::<Vec<_>>().join(",")
            ),
        )]
    };
    (
        StructureMetadata {
            bonds,
            declared_sequences,
        },
        warnings,
    )
}

fn cif_endpoint(
    row: &[String],
    columns: &HashMap<&str, usize>,
    partner: usize,
) -> Option<BondEndpoint> {
    let prefix = format!("_struct_conn.ptnr{partner}_");
    let chain = value(row, columns, &(prefix.clone() + "auth_asym_id"))
        .or_else(|| value(row, columns, &(prefix.clone() + "label_asym_id")))?;
    let residue = value(row, columns, &(prefix.clone() + "auth_seq_id"))
        .or_else(|| value(row, columns, &(prefix.clone() + "label_seq_id")))?;
    let insertion = value(row, columns, &(prefix.clone() + "PDB_ins_code")).unwrap_or("");
    let atom = value(row, columns, &(prefix.clone() + "auth_atom_id"))
        .or_else(|| value(row, columns, &(prefix + "label_atom_id")))?;
    BondEndpoint::new(chain, residue, insertion, atom)
}

fn cif_loops(input: &str) -> ArpeggiaResult<Vec<CifLoop>> {
    let tokens = cif_tokens(input)?;
    let mut loops = Vec::new();
    let mut index = 0;
    while index < tokens.len() {
        if tokens[index] != "loop_" {
            index += 1;
            continue;
        }
        index += 1;
        let mut tags = Vec::new();
        while index < tokens.len() && tokens[index].starts_with('_') {
            tags.push(tokens[index].clone());
            index += 1;
        }
        if tags.is_empty() {
            continue;
        }
        let mut values = Vec::new();
        while index < tokens.len()
            && tokens[index] != "loop_"
            && tokens[index] != "stop_"
            && !tokens[index].starts_with("data_")
            && !tokens[index].starts_with('_')
        {
            values.push(tokens[index].clone());
            index += 1;
        }
        if values.len() % tags.len() != 0 {
            return Err(ArpeggiaError::Parse(format!(
                "mmCIF loop has {} values for {} columns",
                values.len(),
                tags.len()
            )));
        }
        loops.push((
            tags.clone(),
            values.chunks(tags.len()).map(<[String]>::to_vec).collect(),
        ));
    }
    Ok(loops)
}

fn cif_tokens(input: &str) -> ArpeggiaResult<Vec<String>> {
    let mut tokens = Vec::new();
    let mut lines = input.lines();
    while let Some(line) = lines.next() {
        if line.starts_with(';') {
            let mut text = String::new();
            loop {
                let next = lines.next().ok_or_else(|| {
                    ArpeggiaError::Parse("unterminated mmCIF multiline value".into())
                })?;
                if next.starts_with(';') {
                    break;
                }
                if !text.is_empty() {
                    text.push('\n');
                }
                text.push_str(next);
            }
            tokens.push(text);
            continue;
        }
        let mut chars = line.char_indices().peekable();
        while let Some((start, ch)) = chars.next() {
            if ch.is_whitespace() {
                continue;
            }
            if ch == '#' {
                break;
            }
            if ch == '\'' || ch == '"' {
                let quote = ch;
                let value_start = start + ch.len_utf8();
                let mut end = line.len();
                for (position, next) in chars.by_ref() {
                    if next == quote {
                        end = position;
                        break;
                    }
                }
                tokens.push(line[value_start..end].to_string());
            } else {
                let mut end = line.len();
                while let Some(&(position, next)) = chars.peek() {
                    if next.is_whitespace() || next == '#' {
                        end = position;
                        break;
                    }
                    chars.next();
                }
                tokens.push(line[start..end].to_string());
                if chars.peek().is_some_and(|(_, next)| *next == '#') {
                    break;
                }
            }
        }
    }
    Ok(tokens)
}

fn value<'a>(row: &'a [String], columns: &HashMap<&str, usize>, name: &str) -> Option<&'a str> {
    let value = row.get(*columns.get(name)?)?;
    (!matches!(value.as_str(), "." | "?")).then_some(value)
}

fn ordered_bond(first: BondEndpoint, second: BondEndpoint) -> (BondEndpoint, BondEndpoint) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn field(line: &str, start: usize, end: usize) -> &str {
    line.get(start.min(line.len())..end.min(line.len()))
        .unwrap_or("")
}

fn clean(value: &str) -> &str {
    value.trim().trim_matches(['\'', '"'])
}

fn clean_unknown(value: &str) -> &str {
    match clean(value) {
        "." | "?" => "",
        value => value,
    }
}

fn one_letter(name: &str) -> Option<char> {
    Some(match name.trim().to_ascii_uppercase().as_str() {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" | "HID" | "HIE" | "HIP" | "HSD" | "HSE" | "HSP" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" | "MSE" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        "SEC" => 'U',
        "PYL" => 'O',
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pdb_connectivity_and_seqres_are_preserved() {
        let input = "SEQRES   1 A    3  ALA MSE GLY\n\
ATOM      1  SG  CYS A   1       0.000   0.000   0.000  1.00 20.00           S  \n\
ATOM      2  SG  CYS B   2       2.000   0.000   0.000  1.00 20.00           S  \n\
CONECT    1    2\nEND\n";
        let (metadata, warnings) = parse_pdb(input);
        assert!(warnings.is_empty());
        let first = BondEndpoint::new("A", "1", "", "SG").unwrap();
        let second = BondEndpoint::new("B", "2", "", "SG").unwrap();
        assert!(metadata.has_bond(&first, &second));
        assert_eq!(metadata.declared_sequences(), &[("A".into(), "AMG".into())]);
    }

    #[test]
    fn pdb_ssbond_and_link_are_explicit_bonds() {
        let input = "SSBOND   1 CYS A    1    CYS B    2                          1555   1555  2.03\n\
LINK         NZ  LYS A   3                 C1  LIG B   4     1555   1555  1.80\n";
        let (metadata, warnings) = parse_pdb(input);
        assert!(warnings.is_empty());
        assert!(metadata.has_bond(
            &BondEndpoint::from_parts("A", 1, "", "SG"),
            &BondEndpoint::from_parts("B", 2, "", "SG")
        ));
        assert!(metadata.has_bond(
            &BondEndpoint::from_parts("A", 3, "", "NZ"),
            &BondEndpoint::from_parts("B", 4, "", "C1")
        ));
    }

    #[test]
    fn mmcif_shared_entity_maps_to_each_declared_chain() {
        let input = "data_test\n\
loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n1 1 ALA\n1 2 MSE\n\
loop_\n_pdbx_poly_seq_scheme.asym_id\n_pdbx_poly_seq_scheme.entity_id\n_pdbx_poly_seq_scheme.pdb_strand_id\nA 1 H\nB 1 L\n";
        let (metadata, warnings) = parse_mmcif(input).unwrap();
        assert!(warnings.is_empty());
        assert_eq!(
            metadata.declared_sequences(),
            &[("H".into(), "AM".into()), ("L".into(), "AM".into())]
        );
    }

    #[test]
    fn mmcif_struct_conn_covalent_rows_are_explicit_bonds() {
        let input = "data_test\n\
loop_\n_struct_conn.conn_type_id\n_struct_conn.ptnr1_auth_asym_id\n_struct_conn.ptnr1_auth_seq_id\n_struct_conn.ptnr1_auth_atom_id\n_struct_conn.ptnr2_auth_asym_id\n_struct_conn.ptnr2_auth_seq_id\n_struct_conn.ptnr2_auth_atom_id\n\
covale A 1 NZ B 2 C1\nmetalc A 3 ND1 B 4 ZN\n";
        let (metadata, warnings) = parse_mmcif(input).unwrap();
        assert!(warnings.is_empty());
        assert!(metadata.has_bond(
            &BondEndpoint::from_parts("A", 1, "", "NZ"),
            &BondEndpoint::from_parts("B", 2, "", "C1")
        ));
        assert!(!metadata.has_bond(
            &BondEndpoint::from_parts("A", 3, "", "ND1"),
            &BondEndpoint::from_parts("B", 4, "", "ZN")
        ));
    }

    #[test]
    fn declaration_order_and_unknown_monomer_warning_are_preserved() {
        let input = "SEQRES   1 B    1  UNK\nSEQRES   1 A    1  GLY\n";
        let (metadata, warnings) = parse_pdb(input);
        assert_eq!(
            metadata.declared_sequences(),
            &[("B".into(), "X".into()), ("A".into(), "G".into())]
        );
        assert_eq!(warnings.len(), 1);
        assert_eq!(warnings[0].code, WarningCode::UnsupportedMonomer);
    }
}
