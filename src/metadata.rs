//! Narrow bond and declared-sequence parsing missing from `pdbtbx`.

use crate::{Analysis, AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
use pdbtbx::*;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::File;
#[cfg(test)]
use std::io::Cursor;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct ResidueBoundary {
    model: usize,
    chain: String,
    residue: isize,
    insertion: String,
    name: String,
}

impl ResidueBoundary {
    fn from_pdb(model: usize, line: &str) -> Option<Self> {
        Some(Self {
            model,
            chain: clean(field(line, 21, 22)).to_string(),
            residue: clean(field(line, 22, 26)).parse().ok()?,
            insertion: clean_unknown(field(line, 26, 27)).to_string(),
            name: clean(field(line, 17, 20)).to_ascii_uppercase(),
        })
    }
}

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

#[derive(Clone, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
struct BondQualifier {
    model: Option<usize>,
    altloc: Option<String>,
    residue_name: Option<String>,
}

impl BondQualifier {
    fn new(model: Option<usize>, altloc: &str, residue_name: &str) -> Self {
        let optional = |value: &str| {
            let value = clean_unknown(value);
            (!value.is_empty()).then(|| value.to_ascii_uppercase())
        };
        Self {
            model,
            altloc: optional(altloc),
            residue_name: optional(residue_name),
        }
    }

    fn matches(&self, query: &Self) -> bool {
        self.model.is_none_or(|value| query.model == Some(value))
            && self
                .altloc
                .as_ref()
                .is_none_or(|value| query.altloc.as_ref() == Some(value))
            && self
                .residue_name
                .as_ref()
                .is_none_or(|value| query.residue_name.as_ref() == Some(value))
    }
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct QualifiedBondEndpoint {
    endpoint: BondEndpoint,
    qualifier: BondQualifier,
}

impl QualifiedBondEndpoint {
    fn new(endpoint: BondEndpoint, model: Option<usize>, altloc: &str, residue_name: &str) -> Self {
        Self {
            endpoint,
            qualifier: BondQualifier::new(model, altloc, residue_name),
        }
    }

    fn from_entity_atom(entity: &AtomConformerResidueChainModel, atom: &Atom) -> Self {
        Self::new(
            BondEndpoint::from_parts(
                entity.chain().id(),
                entity.residue().serial_number(),
                entity.residue().insertion_code().unwrap_or(""),
                atom.name(),
            ),
            Some(entity.model().serial_number()),
            entity.conformer().alternative_location().unwrap_or(""),
            entity.conformer().name(),
        )
    }
}

type BondKey = (BondEndpoint, BondEndpoint);
type BondRecord = (BondQualifier, BondQualifier, ExplicitBondKind);
type Bonds = HashMap<BondKey, Vec<BondRecord>>;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum ExplicitBondKind {
    Covalent,
    Disulfide,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct AtomIdentity {
    model: usize,
    serial: usize,
}

impl AtomIdentity {
    fn from_entity(entity: &AtomConformerResidueChainModel, atom: &Atom) -> Self {
        Self {
            model: entity.model().serial_number(),
            serial: atom.serial_number(),
        }
    }
}

/// Explicit bonds resolved against the selected atoms of a prepared structure.
#[derive(Clone, Debug, Default)]
pub(crate) struct ResolvedBonds(Vec<(AtomIdentity, AtomIdentity, ExplicitBondKind)>);

impl ResolvedBonds {
    pub(crate) fn contains_entity_atom(
        &self,
        entity: &AtomConformerResidueChainModel,
        atom: &Atom,
    ) -> bool {
        self.kind(
            AtomIdentity::from_entity(entity, entity.atom()),
            AtomIdentity::from_entity(entity, atom),
        )
        .is_some()
    }

    pub(crate) fn kind_entities(
        &self,
        first: &AtomConformerResidueChainModel,
        second: &AtomConformerResidueChainModel,
    ) -> Option<ExplicitBondKind> {
        self.kind(
            AtomIdentity::from_entity(first, first.atom()),
            AtomIdentity::from_entity(second, second.atom()),
        )
    }

    fn kind(&self, first: AtomIdentity, second: AtomIdentity) -> Option<ExplicitBondKind> {
        let atoms = ordered_atom_bond(first, second);
        self.0
            .binary_search_by_key(&atoms, |&(first, second, _)| (first, second))
            .ok()
            .map(|index| self.0[index].2)
    }
}

/// Explicit connectivity and declared polymer sequences from one input file.
#[derive(Clone, Debug, Default)]
pub struct StructureMetadata {
    bonds: Bonds,
    chain_breaks: HashSet<ResidueBoundary>,
    declared_sequences: Vec<(String, String)>,
}

impl StructureMetadata {
    /// Whether the input explicitly connects two atoms.
    pub fn has_bond(&self, first: &BondEndpoint, second: &BondEndpoint) -> bool {
        self.bonds
            .contains_key(&ordered_bond(first.clone(), second.clone()))
    }

    pub(crate) fn resolve_bonds(&self, pdb: &PDB) -> ResolvedBonds {
        if self.bonds.is_empty() {
            return ResolvedBonds::default();
        }
        let mut atoms = HashMap::<BondEndpoint, Vec<(BondQualifier, AtomIdentity)>>::new();
        for entity in pdb.atoms_with_hierarchy() {
            let qualified = QualifiedBondEndpoint::from_entity_atom(&entity, entity.atom());
            atoms.entry(qualified.endpoint).or_default().push((
                qualified.qualifier,
                AtomIdentity::from_entity(&entity, entity.atom()),
            ));
        }

        let mut resolved = Vec::new();
        for ((first, second), records) in &self.bonds {
            let (Some(first_atoms), Some(second_atoms)) = (atoms.get(first), atoms.get(second))
            else {
                continue;
            };
            for (first_record, second_record, kind) in records {
                for (first_qualifier, first_identity) in first_atoms {
                    if !first_record.matches(first_qualifier) {
                        continue;
                    }
                    for (second_qualifier, second_identity) in second_atoms {
                        if second_record.matches(second_qualifier)
                            && first_identity != second_identity
                        {
                            let (first, second) =
                                ordered_atom_bond(*first_identity, *second_identity);
                            resolved.push((first, second, *kind));
                        }
                    }
                }
            }
        }
        resolved.sort_unstable_by_key(|&(first, second, _)| (first, second));
        resolved.dedup_by(|current, previous| {
            let same_atoms = current.0 == previous.0 && current.1 == previous.1;
            if same_atoms
                && (current.2 == ExplicitBondKind::Disulfide
                    || previous.2 == ExplicitBondKind::Disulfide)
            {
                previous.2 = ExplicitBondKind::Disulfide;
            }
            same_atoms
        });
        ResolvedBonds(resolved)
    }

    #[cfg(test)]
    fn has_qualified_bond(
        &self,
        first: QualifiedBondEndpoint,
        second: QualifiedBondEndpoint,
    ) -> bool {
        let (first, second) = ordered_qualified_bond(first, second);
        self.bonds
            .get(&(first.endpoint, second.endpoint))
            .is_some_and(|records| {
                records.iter().any(|(record_first, record_second, _)| {
                    record_first.matches(&first.qualifier)
                        && record_second.matches(&second.qualifier)
                })
            })
    }

    /// Chain-oriented declared polymer sequences.
    pub fn declared_sequences(&self) -> &[(String, String)] {
        &self.declared_sequences
    }

    pub(crate) fn has_chain_break_after(
        &self,
        model: usize,
        chain: &str,
        residue: isize,
        insertion: &str,
        name: &str,
    ) -> bool {
        self.chain_breaks.contains(&ResidueBoundary {
            model,
            chain: chain.to_string(),
            residue,
            insertion: insertion.to_string(),
            name: name.to_ascii_uppercase(),
        })
    }
}

/// Read only connectivity and polymer-sequence metadata from PDB or mmCIF.
pub fn read_metadata(path: impl AsRef<Path>) -> ArpeggiaResult<Analysis<StructureMetadata>> {
    let path = path.as_ref();
    let (metadata, warnings) = if is_mmcif(path)? {
        parse_mmcif_reader(BufReader::new(File::open(path)?))?
    } else {
        let input = std::fs::read_to_string(path)?;
        parse_pdb(&input)?
    };
    Ok(Analysis::new(metadata, warnings))
}

/// Return declared PDB `SEQRES` or mmCIF entity-polymer sequences by chain.
pub fn get_seqres(path: impl AsRef<Path>) -> ArpeggiaResult<Analysis<Vec<(String, String)>>> {
    let path = path.as_ref();
    if is_mmcif(path)? {
        let (metadata, warnings) = parse_mmcif_reader(BufReader::new(File::open(path)?))?;
        return Ok(Analysis::new(metadata.declared_sequences, warnings));
    }
    let input = std::fs::read_to_string(path)?;
    let sequences = pdb_declared_sequences(&input);
    let (metadata, warnings) = metadata_with_sequences(HashMap::new(), HashSet::new(), sequences);
    Ok(Analysis::new(metadata.declared_sequences, warnings))
}

fn is_mmcif(path: &Path) -> ArpeggiaResult<bool> {
    if path
        .extension()
        .and_then(|value| value.to_str())
        .is_some_and(|extension| matches!(extension.to_ascii_lowercase().as_str(), "cif" | "mmcif"))
    {
        return Ok(true);
    }
    let mut reader = BufReader::new(File::open(path)?);
    let prefix = reader.fill_buf()?;
    Ok(std::str::from_utf8(prefix)
        .unwrap_or("")
        .trim_start()
        .starts_with("data_"))
}

fn parse_pdb(input: &str) -> ArpeggiaResult<(StructureMetadata, Vec<AnalysisWarning>)> {
    let mut bonds = HashMap::new();
    let mut chain_breaks = HashSet::new();
    let mut atoms = HashMap::<usize, Vec<QualifiedBondEndpoint>>::new();
    let mut conect = Vec::new();
    let mut current_model = 0;
    let mut last_coordinate = None;

    for line in input.lines() {
        match field(line, 0, 6).trim() {
            "MODEL" => {
                current_model = clean(field(line, 10, 14)).parse().unwrap_or(current_model);
                last_coordinate = None;
            }
            "ATOM" | "HETATM" => {
                let residue = ResidueBoundary::from_pdb(current_model, line).ok_or_else(|| {
                    ArpeggiaError::Parse("invalid residue identity in coordinate record".into())
                })?;
                if chain_breaks.contains(&residue) {
                    return Err(ArpeggiaError::Parse(format!(
                        "TER occurs inside model {}, chain {}, residue {}{}",
                        residue.model, residue.chain, residue.residue, residue.insertion
                    )));
                }
                last_coordinate = Some(residue);
                if let (Ok(serial), Some(endpoint)) = (
                    field(line, 6, 11).trim().parse::<usize>(),
                    BondEndpoint::new(
                        field(line, 21, 22),
                        field(line, 22, 26),
                        field(line, 26, 27),
                        field(line, 12, 16),
                    ),
                ) {
                    atoms
                        .entry(serial)
                        .or_default()
                        .push(QualifiedBondEndpoint::new(
                            endpoint,
                            Some(current_model),
                            field(line, 16, 17),
                            field(line, 17, 20),
                        ));
                }
            }
            "TER" => {
                let residue = ResidueBoundary::from_pdb(current_model, line)
                    .ok_or_else(|| ArpeggiaError::Parse("invalid TER residue identity".into()))?;
                if last_coordinate.as_ref() != Some(&residue) {
                    return Err(ArpeggiaError::Parse(format!(
                        "TER does not identify the preceding complete residue in model {}, chain {}",
                        current_model, residue.chain
                    )));
                }
                chain_breaks.insert(residue);
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
                    insert_bond(
                        &mut bonds,
                        QualifiedBondEndpoint::new(first, None, "", "CYS"),
                        QualifiedBondEndpoint::new(second, None, "", "CYS"),
                        ExplicitBondKind::Disulfide,
                    );
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
                    insert_bond(
                        &mut bonds,
                        QualifiedBondEndpoint::new(
                            first,
                            None,
                            field(line, 16, 17),
                            field(line, 17, 20),
                        ),
                        QualifiedBondEndpoint::new(
                            second,
                            None,
                            field(line, 46, 47),
                            field(line, 47, 50),
                        ),
                        ExplicitBondKind::Covalent,
                    );
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
            _ => {}
        }
    }
    for (first, second) in conect {
        if let (Some(first), Some(second)) = (atoms.get(&first), atoms.get(&second)) {
            for first in first {
                for second in second {
                    if first.qualifier.model == second.qualifier.model {
                        insert_bond(
                            &mut bonds,
                            first.clone(),
                            second.clone(),
                            ExplicitBondKind::Covalent,
                        );
                    }
                }
            }
        }
    }
    Ok(metadata_with_sequences(
        bonds,
        chain_breaks,
        pdb_declared_sequences(input),
    ))
}

fn pdb_declared_sequences(input: &str) -> Vec<(String, Vec<String>)> {
    let mut sequences = Vec::<(String, Vec<String>)>::new();
    for line in input
        .lines()
        .filter(|line| field(line, 0, 6).trim() == "SEQRES")
    {
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
    sequences
}

#[cfg(test)]
fn parse_mmcif(input: &str) -> ArpeggiaResult<(StructureMetadata, Vec<AnalysisWarning>)> {
    parse_mmcif_reader(Cursor::new(input))
}

fn parse_mmcif_reader(
    reader: impl BufRead,
) -> ArpeggiaResult<(StructureMetadata, Vec<AnalysisWarning>)> {
    let mut tokens = CifTokens::new(reader);
    let mut pending = None;
    let mut bonds = HashMap::new();
    let mut entities = BTreeMap::<String, BTreeMap<usize, String>>::new();
    let mut chain_entities = Vec::<(String, String)>::new();

    while let Some(token) = next_cif_token(&mut tokens, &mut pending, true)? {
        if !token.is_keyword("loop_") {
            continue;
        }
        let mut tags = Vec::new();
        while let Some(token) = tokens.next(true)? {
            if token.bare && token.value.starts_with('_') {
                tags.push(token.value);
            } else {
                pending = Some(token);
                break;
            }
        }
        if tags.is_empty() {
            continue;
        }
        let columns: HashMap<&str, usize> = tags
            .iter()
            .enumerate()
            .map(|(index, tag)| (tag.as_str(), index))
            .collect();
        let target = if tags.iter().any(|tag| tag.starts_with("_struct_conn.")) {
            CifTarget::StructConn
        } else if tags.iter().any(|tag| tag.starts_with("_entity_poly_seq.")) {
            CifTarget::EntityPolySeq
        } else if tags
            .iter()
            .any(|tag| tag.starts_with("_pdbx_poly_seq_scheme."))
        {
            CifTarget::PolySeqScheme
        } else {
            CifTarget::Ignore
        };
        let mut row = Vec::with_capacity(if target == CifTarget::Ignore {
            0
        } else {
            tags.len()
        });
        let mut value_count = 0;

        while let Some(token) =
            next_cif_token(&mut tokens, &mut pending, target != CifTarget::Ignore)?
        {
            if token.is_control() {
                pending = Some(token);
                break;
            }
            value_count += 1;
            if target != CifTarget::Ignore {
                row.push(token.value);
            }
            if value_count % tags.len() != 0 {
                continue;
            }

            match target {
                CifTarget::StructConn => {
                    let kind = value(&row, &columns, "_struct_conn.conn_type_id").unwrap_or("");
                    if (kind.starts_with("covale") || kind.starts_with("disulf"))
                        && let (Some(first), Some(second)) = (
                            cif_endpoint(&row, &columns, 1),
                            cif_endpoint(&row, &columns, 2),
                        )
                    {
                        let kind = if kind.starts_with("disulf") {
                            ExplicitBondKind::Disulfide
                        } else {
                            ExplicitBondKind::Covalent
                        };
                        insert_bond(&mut bonds, first, second, kind);
                    }
                }
                CifTarget::EntityPolySeq => {
                    if let (Some(entity), Some(number), Some(monomer)) = (
                        value(&row, &columns, "_entity_poly_seq.entity_id"),
                        value(&row, &columns, "_entity_poly_seq.num")
                            .and_then(|number| number.parse::<usize>().ok()),
                        value(&row, &columns, "_entity_poly_seq.mon_id"),
                    ) {
                        entities
                            .entry(entity.to_string())
                            .or_default()
                            .insert(number, monomer.to_string());
                    }
                }
                CifTarget::PolySeqScheme => {
                    if let (Some(chain), Some(entity)) = (
                        value(&row, &columns, "_pdbx_poly_seq_scheme.pdb_strand_id")
                            .or_else(|| value(&row, &columns, "_pdbx_poly_seq_scheme.asym_id")),
                        value(&row, &columns, "_pdbx_poly_seq_scheme.entity_id"),
                    ) {
                        for chain in chain.split(',') {
                            let chain = clean(chain).to_string();
                            if !chain_entities.iter().any(|(known, _)| known == &chain) {
                                chain_entities.push((chain, entity.to_string()));
                            }
                        }
                    }
                }
                CifTarget::Ignore => {}
            }

            row.clear();
        }
        if value_count % tags.len() != 0 {
            return Err(ArpeggiaError::Parse(format!(
                "mmCIF loop has {} values for {} columns",
                value_count,
                tags.len()
            )));
        }
    }

    let sequences = chain_entities
        .into_iter()
        .filter_map(|(chain, entity)| {
            let monomers = entities.get(&entity)?;
            Some((chain, monomers.values().cloned().collect::<Vec<_>>()))
        })
        .collect();
    Ok(metadata_with_sequences(bonds, HashSet::new(), sequences))
}

fn metadata_with_sequences(
    bonds: Bonds,
    chain_breaks: HashSet<ResidueBoundary>,
    sequences: Vec<(String, Vec<String>)>,
) -> (StructureMetadata, Vec<AnalysisWarning>) {
    let mut unsupported = BTreeSet::new();
    let declared_sequences = sequences
        .into_iter()
        .map(|(chain, monomers)| {
            let sequence = monomers
                .iter()
                .map(
                    |name| match crate::contacts::residues::one_letter_code(name) {
                        Some(letter) => letter.chars().next().unwrap(),
                        None => {
                            unsupported.insert(name.trim().to_ascii_uppercase());
                            'X'
                        }
                    },
                )
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
            chain_breaks,
            declared_sequences,
        },
        warnings,
    )
}

fn cif_endpoint(
    row: &[String],
    columns: &HashMap<&str, usize>,
    partner: usize,
) -> Option<QualifiedBondEndpoint> {
    let prefix = format!("_struct_conn.ptnr{partner}_");
    let pdbx_prefix = format!("_struct_conn.pdbx_ptnr{partner}_");
    let chain = value(row, columns, &(prefix.clone() + "auth_asym_id"))
        .or_else(|| value(row, columns, &(prefix.clone() + "label_asym_id")))?;
    let residue = value(row, columns, &(prefix.clone() + "auth_seq_id"))
        .or_else(|| value(row, columns, &(prefix.clone() + "label_seq_id")))?;
    let insertion = value(row, columns, &(pdbx_prefix.clone() + "PDB_ins_code")).unwrap_or("");
    let atom = value(row, columns, &(prefix.clone() + "label_atom_id"))
        .or_else(|| value(row, columns, &(prefix.clone() + "auth_atom_id")))?;
    let residue_name = value(row, columns, &(prefix.clone() + "label_comp_id"))
        .or_else(|| value(row, columns, &(prefix + "auth_comp_id")))
        .unwrap_or("");
    let altloc = value(row, columns, &(pdbx_prefix.clone() + "label_alt_id")).unwrap_or("");
    let model =
        value(row, columns, &(pdbx_prefix + "PDB_model_num")).and_then(|value| value.parse().ok());
    Some(QualifiedBondEndpoint::new(
        BondEndpoint::new(chain, residue, insertion, atom)?,
        model,
        altloc,
        residue_name,
    ))
}

#[derive(Clone, Copy, Eq, PartialEq)]
enum CifTarget {
    Ignore,
    StructConn,
    EntityPolySeq,
    PolySeqScheme,
}

#[derive(Debug)]
struct CifToken {
    value: String,
    bare: bool,
}

impl CifToken {
    fn is_keyword(&self, keyword: &str) -> bool {
        self.bare && self.value.eq_ignore_ascii_case(keyword)
    }

    fn is_control(&self) -> bool {
        self.bare && is_cif_control(&self.value)
    }
}

fn is_cif_control(value: &str) -> bool {
    value.starts_with('_')
        || value.eq_ignore_ascii_case("loop_")
        || value.eq_ignore_ascii_case("stop_")
        || has_ascii_prefix(value, "data_")
        || has_ascii_prefix(value, "save_")
        || value.eq_ignore_ascii_case("global_")
}

fn has_ascii_prefix(value: &str, prefix: &str) -> bool {
    value
        .get(..prefix.len())
        .is_some_and(|candidate| candidate.eq_ignore_ascii_case(prefix))
}

struct CifTokens<R> {
    reader: R,
    line: String,
    offset: usize,
}

impl<R: BufRead> CifTokens<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            line: String::new(),
            offset: 0,
        }
    }

    fn next(&mut self, capture: bool) -> ArpeggiaResult<Option<CifToken>> {
        loop {
            if let Some(token) = self.line_token(capture)? {
                return Ok(Some(token));
            }
            self.line.clear();
            self.offset = 0;
            if self.reader.read_line(&mut self.line)? == 0 {
                return Ok(None);
            }
            if self.line.starts_with(';') {
                return self.multiline_value(capture).map(Some);
            }
        }
    }

    fn line_token(&mut self, capture: bool) -> ArpeggiaResult<Option<CifToken>> {
        let mut chars = self.line[self.offset..].char_indices().peekable();
        while let Some((relative_start, ch)) = chars.next() {
            let start = self.offset + relative_start;
            if ch.is_whitespace() {
                continue;
            }
            if ch == '#' {
                self.offset = self.line.len();
                return Ok(None);
            }
            if ch == '\'' || ch == '"' {
                let quote = ch;
                let value_start = start + ch.len_utf8();
                for (relative_end, next) in chars.by_ref() {
                    let end = self.offset + relative_end;
                    let after_quote = end + next.len_utf8();
                    if next == quote
                        && self.line[after_quote..]
                            .chars()
                            .next()
                            .is_none_or(|next| next.is_whitespace() || next == '#')
                    {
                        self.offset = after_quote;
                        return Ok(Some(CifToken {
                            value: if capture {
                                self.line[value_start..end].to_string()
                            } else {
                                String::new()
                            },
                            bare: false,
                        }));
                    }
                }
                return Err(ArpeggiaError::Parse(
                    "unterminated quoted mmCIF value".into(),
                ));
            }

            let mut end = self.line.len();
            while let Some(&(relative_end, next)) = chars.peek() {
                if next.is_whitespace() {
                    end = self.offset + relative_end;
                    break;
                }
                chars.next();
            }
            self.offset = end;
            let value = self.line[start..end].trim_end_matches(['\r', '\n']);
            let control = is_cif_control(value);
            return Ok(Some(CifToken {
                value: if capture || control {
                    value.to_string()
                } else {
                    String::new()
                },
                bare: true,
            }));
        }
        self.offset = self.line.len();
        Ok(None)
    }

    fn multiline_value(&mut self, capture: bool) -> ArpeggiaResult<CifToken> {
        let mut value = if capture {
            self.line[1..].trim_end_matches(['\r', '\n']).to_string()
        } else {
            String::new()
        };
        loop {
            self.line.clear();
            if self.reader.read_line(&mut self.line)? == 0 {
                return Err(ArpeggiaError::Parse(
                    "unterminated mmCIF multiline value".into(),
                ));
            }
            if self.line.starts_with(';') {
                self.offset = self.line.len();
                return Ok(CifToken { value, bare: false });
            }
            if capture {
                if !value.is_empty() {
                    value.push('\n');
                }
                value.push_str(self.line.trim_end_matches(['\r', '\n']));
            }
        }
    }
}

fn next_cif_token<R: BufRead>(
    tokens: &mut CifTokens<R>,
    pending: &mut Option<CifToken>,
    capture: bool,
) -> ArpeggiaResult<Option<CifToken>> {
    if pending.is_some() {
        Ok(pending.take())
    } else {
        tokens.next(capture)
    }
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

fn ordered_atom_bond(first: AtomIdentity, second: AtomIdentity) -> (AtomIdentity, AtomIdentity) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn ordered_qualified_bond(
    first: QualifiedBondEndpoint,
    second: QualifiedBondEndpoint,
) -> (QualifiedBondEndpoint, QualifiedBondEndpoint) {
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

fn insert_bond(
    bonds: &mut Bonds,
    first: QualifiedBondEndpoint,
    second: QualifiedBondEndpoint,
    kind: ExplicitBondKind,
) {
    let (first, second) = ordered_qualified_bond(first, second);
    let record = (first.qualifier, second.qualifier, kind);
    let records = bonds.entry((first.endpoint, second.endpoint)).or_default();
    if !records.contains(&record) {
        records.push(record);
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pdb_connectivity_and_seqres_are_preserved() {
        let input = "SEQRES   1 A    3  ALA MSE GLY\n\
ATOM      1  SG  CYS A   1       0.000   0.000   0.000  1.00 20.00           S  \n\
ATOM      2  SG  CYS B   2       2.000   0.000   0.000  1.00 20.00           S  \n\
CONECT    1    2\nEND\n";
        let (metadata, warnings) = parse_pdb(input).unwrap();
        assert!(warnings.is_empty());
        let first = BondEndpoint::new("A", "1", "", "SG").unwrap();
        let second = BondEndpoint::new("B", "2", "", "SG").unwrap();
        assert!(metadata.has_bond(&first, &second));
        assert_eq!(metadata.declared_sequences(), &[("A".into(), "AMG".into())]);
    }

    #[test]
    fn malformed_ter_fails_metadata_but_not_seqres() {
        let input = "SEQRES   1 A    1  ALA\n\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N  \n\
TER       2      ALA A   1\n\
ATOM      3  CA  ALA A   1       1.400   0.000   0.000  1.00 20.00           C  \n\
END\n";
        let path =
            std::env::temp_dir().join(format!("arpeggia-malformed-ter-{}.pdb", std::process::id()));
        std::fs::write(&path, input).unwrap();

        let metadata = read_metadata(&path);
        let seqres = get_seqres(&path).unwrap();
        std::fs::remove_file(path).unwrap();

        assert!(matches!(metadata, Err(ArpeggiaError::Parse(_))));
        assert_eq!(seqres.value, vec![("A".into(), "A".into())]);
    }

    #[test]
    fn pdb_ssbond_and_link_are_explicit_bonds() {
        let input = "SSBOND   1 CYS A    1    CYS B    2                          1555   1555  2.03\n\
LINK         NZ  LYS A   3                 C1  LIG B   4     1555   1555  1.80\n";
        let (metadata, warnings) = parse_pdb(input).unwrap();
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
    fn pdb_link_applies_only_to_its_named_altloc() {
        let input =
            "LINK         SG BCYS A   1                 SG BCYS B   2     1555   1555  2.03\n";
        let (metadata, _) = parse_pdb(input).unwrap();
        let endpoint = |chain, residue, altloc| {
            QualifiedBondEndpoint::new(
                BondEndpoint::from_parts(chain, residue, "", "SG"),
                Some(1),
                altloc,
                "CYS",
            )
        };

        assert!(metadata.has_qualified_bond(endpoint("A", 1, "B"), endpoint("B", 2, "B")));
        assert!(!metadata.has_qualified_bond(endpoint("A", 1, "A"), endpoint("B", 2, "A")));
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
loop_\n_struct_conn.conn_type_id\n\
_struct_conn.ptnr1_auth_asym_id\n_struct_conn.ptnr1_auth_seq_id\n\
_struct_conn.ptnr1_label_atom_id\n_struct_conn.ptnr1_auth_atom_id\n\
_struct_conn.ptnr1_label_comp_id\n_struct_conn.pdbx_ptnr1_label_alt_id\n\
_struct_conn.pdbx_ptnr1_PDB_ins_code\n\
_struct_conn.ptnr2_auth_asym_id\n_struct_conn.ptnr2_auth_seq_id\n\
_struct_conn.ptnr2_label_atom_id\n_struct_conn.ptnr2_auth_atom_id\n\
_struct_conn.ptnr2_label_comp_id\n_struct_conn.pdbx_ptnr2_label_alt_id\n\
_struct_conn.pdbx_ptnr2_PDB_ins_code\n\
covale A 1 NZ NZ_AUTH LYS A I B 2 C1 C1_AUTH LIG B ?\n\
metalc A 3 ND1 ND1_AUTH HIS . ? B 4 ZN ZN_AUTH ZN . ?\n";
        let (metadata, warnings) = parse_mmcif(input).unwrap();
        assert!(warnings.is_empty());
        assert!(metadata.has_bond(
            &BondEndpoint::from_parts("A", 1, "I", "NZ"),
            &BondEndpoint::from_parts("B", 2, "", "C1")
        ));
        assert!(!metadata.has_bond(
            &BondEndpoint::from_parts("A", 3, "", "ND1"),
            &BondEndpoint::from_parts("B", 4, "", "ZN")
        ));
        assert!(metadata.has_qualified_bond(
            QualifiedBondEndpoint::new(
                BondEndpoint::from_parts("A", 1, "I", "NZ"),
                Some(1),
                "A",
                "LYS"
            ),
            QualifiedBondEndpoint::new(
                BondEndpoint::from_parts("B", 2, "", "C1"),
                Some(1),
                "B",
                "LIG"
            )
        ));
        assert!(!metadata.has_qualified_bond(
            QualifiedBondEndpoint::new(
                BondEndpoint::from_parts("A", 1, "I", "NZ"),
                Some(1),
                "B",
                "LYS"
            ),
            QualifiedBondEndpoint::new(
                BondEndpoint::from_parts("B", 2, "", "C1"),
                Some(1),
                "B",
                "LIG"
            )
        ));
    }

    #[test]
    fn mmcif_quoted_reserved_value_is_not_loop_syntax() {
        let input = "data_test\n\
loop_\n_entity_poly_seq.entity_id\n_entity_poly_seq.num\n_entity_poly_seq.mon_id\n\
1 1 '_UNK'\n\
loop_\n_pdbx_poly_seq_scheme.asym_id\n_pdbx_poly_seq_scheme.entity_id\n\
A 1\n";
        let path = std::env::temp_dir().join(format!(
            "arpeggia-quoted-cif-value-{}.cif",
            std::process::id()
        ));
        std::fs::write(&path, input).unwrap();

        let sequences = get_seqres(&path);
        std::fs::remove_file(path).unwrap();

        let analysis = sequences.unwrap();
        assert_eq!(analysis.value, vec![("A".into(), "X".into())]);
        assert_eq!(analysis.warnings[0].code, WarningCode::UnsupportedMonomer);
    }

    #[test]
    fn declaration_order_and_unknown_monomer_warning_are_preserved() {
        let input = "SEQRES   1 B    1  UNK\nSEQRES   1 A    1  GLY\n";
        let (metadata, warnings) = parse_pdb(input).unwrap();
        assert_eq!(
            metadata.declared_sequences(),
            &[("B".into(), "X".into()), ("A".into(), "G".into())]
        );
        assert_eq!(warnings.len(), 1);
        assert_eq!(warnings[0].code, WarningCode::UnsupportedMonomer);
    }
}
