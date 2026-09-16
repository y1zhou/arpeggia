//! Residue selectors, atom populations, and exact coordinate correspondence.

use crate::contacts::one_letter_code;
use crate::structure::selected_model;
use crate::{AnalysisWarning, ArpeggiaError, ArpeggiaResult, WarningCode};
use pdbtbx::{ContainsAtomConformer, Element, PDB};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

/// Atom population used for structural superposition.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, clap::ValueEnum, Serialize)]
pub enum AtomSubset {
    /// Alpha-carbon atoms only.
    #[default]
    Ca,
    /// Protein N, CA, C, O, and OXT atoms.
    Backbone,
    /// All non-hydrogen atoms in protein residues and terminal caps.
    Heavy,
    /// All atoms in protein residues and terminal caps.
    All,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(super) struct AtomIdentity {
    pub(super) chain: String,
    pub(super) residue_number: isize,
    pub(super) insertion: String,
    pub(super) residue_name: String,
    pub(super) atom_name: String,
}

#[derive(Clone, Debug)]
pub(super) struct SelectedCoordinateUnion {
    pub(super) keys: Vec<AtomIdentity>,
    pub(super) coordinates: Vec<[f64; 3]>,
    pub(super) superpose_end: usize,
    pub(super) rmsd_start: usize,
    pub(super) warnings: Vec<AnalysisWarning>,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub(super) enum SelectionMembership {
    SuperposeOnly,
    Both,
    RmsdOnly,
}

pub(super) fn select_coordinate_union(
    pdb: &PDB,
    model_num: usize,
    superpose_selector: &ResidueSelector,
    rmsd_selector: &ResidueSelector,
    atoms: AtomSubset,
) -> ArpeggiaResult<SelectedCoordinateUnion> {
    let model = selected_model(pdb, model_num)?;
    superpose_selector
        .validate_chains(model.chains().map(|chain| chain.id()))
        .map_err(|error| selection_error("Superposition Selection", error))?;
    rmsd_selector
        .validate_chains(model.chains().map(|chain| chain.id()))
        .map_err(|error| selection_error("RMSD Selection", error))?;

    let mut selected = Vec::new();
    for chain in model.chains() {
        for residue in chain.residues() {
            let superpose = superpose_selector.matches(
                chain.id(),
                residue.serial_number(),
                residue.insertion_code(),
            );
            let rmsd = rmsd_selector.matches(
                chain.id(),
                residue.serial_number(),
                residue.insertion_code(),
            );
            let membership = match (superpose, rmsd) {
                (true, true) => SelectionMembership::Both,
                (true, false) => SelectionMembership::SuperposeOnly,
                (false, true) => SelectionMembership::RmsdOnly,
                (false, false) => continue,
            };
            let residue_name = residue.name().unwrap_or("");
            // TODO: Add an exact-correspondence all-atom mode that retains arbitrary
            // polymers, ligands, and modified residues when both inputs match.
            let is_protein = residue.name().and_then(one_letter_code).is_some();
            let is_cap = matches!(residue_name, "ACE" | "NH2");
            if !(is_protein || is_cap && matches!(atoms, AtomSubset::Heavy | AtomSubset::All)) {
                continue;
            }
            for hierarchy in residue.atoms_with_hierarchy() {
                let atom = hierarchy.atom();
                if !atom_in_subset(atom.name(), atom.element(), atoms) {
                    continue;
                }
                let coordinate = [atom.x(), atom.y(), atom.z()];
                if !coordinate.iter().all(|value| value.is_finite()) {
                    return Err(ArpeggiaError::Calculation(format!(
                        "atom {} {}:{}{} {} has non-finite coordinates",
                        atom.name(),
                        chain.id(),
                        residue.serial_number(),
                        residue.insertion_code().unwrap_or(""),
                        residue_name
                    )));
                }
                selected.push((
                    membership,
                    AtomIdentity {
                        chain: chain.id().to_string(),
                        residue_number: residue.serial_number(),
                        insertion: residue.insertion_code().unwrap_or("").to_string(),
                        residue_name: residue_name.to_string(),
                        atom_name: atom.name().to_string(),
                    },
                    coordinate,
                ));
            }
        }
    }
    selected.sort_unstable_by(|left, right| (&left.0, &left.1).cmp(&(&right.0, &right.1)));
    if let Some(duplicate) = selected
        .windows(2)
        .find(|pair| pair[0].1 == pair[1].1)
        .map(|pair| &pair[0].1)
    {
        return Err(ArpeggiaError::Calculation(format!(
            "selected atom identity occurs more than once: {duplicate:?}"
        )));
    }
    let superpose_end =
        selected.partition_point(|(membership, _, _)| *membership != SelectionMembership::RmsdOnly);
    let rmsd_start = selected
        .partition_point(|(membership, _, _)| *membership == SelectionMembership::SuperposeOnly);
    if superpose_end == 0 {
        return Err(ArpeggiaError::InvalidArgument(
            "Superposition Selection contains no atoms".into(),
        ));
    }
    if rmsd_start == selected.len() {
        return Err(ArpeggiaError::InvalidArgument(
            "RMSD Selection contains no atoms".into(),
        ));
    }

    let mut warnings = Vec::new();
    if model_num == 0 && pdb.model_count() > 1 {
        warnings.push(AnalysisWarning::new(
            WarningCode::ModelSelected,
            format!(
                "selected first model {} from {} coordinate models",
                model.serial_number(),
                pdb.model_count()
            ),
        ));
    }
    let (keys, coordinates) = selected
        .into_iter()
        .map(|(_, identity, coordinate)| (identity, coordinate))
        .unzip();
    Ok(SelectedCoordinateUnion {
        keys,
        coordinates,
        superpose_end,
        rmsd_start,
        warnings,
    })
}

fn selection_error(selection: &str, error: ArpeggiaError) -> ArpeggiaError {
    match error {
        ArpeggiaError::InvalidArgument(message) => {
            ArpeggiaError::InvalidArgument(format!("{selection}: {message}"))
        }
        error => error,
    }
}

pub(super) fn atom_in_subset(name: &str, element: Option<&Element>, subset: AtomSubset) -> bool {
    match subset {
        AtomSubset::Ca => name == "CA",
        AtomSubset::Backbone => matches!(name, "N" | "CA" | "C" | "O" | "OXT"),
        AtomSubset::Heavy => {
            element != Some(&Element::H) && (element.is_some() || !is_hydrogen_atom_name(name))
        }
        AtomSubset::All => true,
    }
}

fn is_hydrogen_atom_name(name: &str) -> bool {
    matches!(
        name.trim_start_matches(|character: char| character.is_ascii_digit())
            .as_bytes()
            .first(),
        Some(b'H' | b'D' | b'h' | b'd')
    )
}

pub(super) fn validate_correspondence(
    reference: &[AtomIdentity],
    query: &[AtomIdentity],
) -> ArpeggiaResult<()> {
    if reference.len() != query.len() {
        return Err(ArpeggiaError::Calculation(format!(
            "atom correspondence mismatch: reference has {} selected atoms but query has {}",
            reference.len(),
            query.len()
        )));
    }
    if let Some((index, (left, right))) = reference
        .iter()
        .zip(query)
        .enumerate()
        .find(|(_, (left, right))| left != right)
    {
        return Err(ArpeggiaError::Calculation(format!(
            "atom correspondence mismatch at position {index}: {left:?} != {right:?}"
        )));
    }
    Ok(())
}

pub(super) fn validate_selection_correspondence(
    reference: &SelectedCoordinateUnion,
    query: &SelectedCoordinateUnion,
) -> ArpeggiaResult<()> {
    validate_selection_keys(
        &reference.keys,
        reference.superpose_end,
        reference.rmsd_start,
        &query.keys,
        query.superpose_end,
        query.rmsd_start,
    )
}

pub(super) fn validate_selection_keys(
    reference: &[AtomIdentity],
    reference_superpose_end: usize,
    reference_rmsd_start: usize,
    query: &[AtomIdentity],
    query_superpose_end: usize,
    query_rmsd_start: usize,
) -> ArpeggiaResult<()> {
    validate_correspondence(
        &reference[..reference_superpose_end],
        &query[..query_superpose_end],
    )
    .map_err(|error| selection_error("Superposition Selection", error))?;
    validate_correspondence(
        &reference[reference_rmsd_start..],
        &query[query_rmsd_start..],
    )
    .map_err(|error| selection_error("RMSD Selection", error))
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct ResidueBound {
    number: isize,
    insertion: Option<String>,
}

impl ResidueBound {
    fn parse(value: &str) -> Option<Self> {
        let digit_start = usize::from(value.starts_with('-'));
        let digit_count = value[digit_start..]
            .bytes()
            .take_while(u8::is_ascii_digit)
            .count();
        if digit_count == 0 {
            return None;
        }
        let number_end = digit_start + digit_count;
        let number = value[..number_end].parse().ok()?;
        let insertion = &value[number_end..];
        if insertion.is_empty() {
            Some(Self {
                number,
                insertion: None,
            })
        } else if insertion.len() == 1 && insertion.bytes().all(|byte| byte.is_ascii_alphanumeric())
        {
            Some(Self {
                number,
                insertion: Some(insertion.to_ascii_uppercase()),
            })
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct ResidueSpan {
    start: ResidueBound,
    end: ResidueBound,
}

impl ResidueSpan {
    fn parse(value: &str) -> ArpeggiaResult<Self> {
        if let Some(bound) = ResidueBound::parse(value) {
            return Ok(Self {
                start: bound.clone(),
                end: bound,
            });
        }
        for (index, byte) in value.bytes().enumerate().skip(1) {
            if byte != b'-' {
                continue;
            }
            let Some(start) = ResidueBound::parse(&value[..index]) else {
                continue;
            };
            let Some(end) = ResidueBound::parse(&value[index + 1..]) else {
                continue;
            };
            if start.number > end.number
                || (start.number == end.number
                    && end.insertion.as_ref().is_some_and(|end| {
                        start.insertion.as_ref().is_some_and(|start| start > end)
                    }))
            {
                return Err(ArpeggiaError::InvalidArgument(format!(
                    "residue range starts after it ends: {value}"
                )));
            }
            return Ok(Self { start, end });
        }
        Err(ArpeggiaError::InvalidArgument(format!(
            "invalid residue selection: {value}"
        )))
    }

    fn matches(&self, number: isize, insertion: Option<&str>) -> bool {
        lower_bound_matches(number, insertion, &self.start)
            && upper_bound_matches(number, insertion, &self.end)
    }
}

fn lower_bound_matches(number: isize, insertion: Option<&str>, bound: &ResidueBound) -> bool {
    number > bound.number
        || (number == bound.number
            && bound
                .insertion
                .as_deref()
                .is_none_or(|start| insertion.unwrap_or("") >= start))
}

fn upper_bound_matches(number: isize, insertion: Option<&str>, bound: &ResidueBound) -> bool {
    number < bound.number
        || (number == bound.number
            && bound
                .insertion
                .as_deref()
                .is_none_or(|end| insertion.unwrap_or("") <= end))
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub(super) struct ResidueSelector {
    chains: BTreeMap<String, Option<Vec<ResidueSpan>>>,
}

impl ResidueSelector {
    pub(super) fn parse(value: &str) -> ArpeggiaResult<Self> {
        if value.trim().is_empty() {
            return Ok(Self::default());
        }
        let mut chains: BTreeMap<String, Option<Vec<ResidueSpan>>> = BTreeMap::new();
        for clause in value.split(',').map(str::trim) {
            if clause.is_empty() {
                return Err(ArpeggiaError::InvalidArgument(
                    "residue selection contains an empty clause".into(),
                ));
            }
            let (chain, span) = clause
                .split_once(':')
                .map_or((clause, None), |(chain, span)| (chain, Some(span)));
            if chain.is_empty() {
                return Err(ArpeggiaError::InvalidArgument(
                    "residue selection contains an empty chain".into(),
                ));
            }
            match span {
                None => {
                    chains.insert(chain.to_string(), None);
                }
                Some("") => {
                    return Err(ArpeggiaError::InvalidArgument(format!(
                        "chain {chain} has an empty residue selection"
                    )));
                }
                Some(span) => {
                    if chains.get(chain).is_some_and(Option::is_none) {
                        continue;
                    }
                    chains
                        .entry(chain.to_string())
                        .or_insert_with(|| Some(Vec::new()))
                        .as_mut()
                        .expect("entry was initialized with a residue list")
                        .push(ResidueSpan::parse(span)?);
                }
            }
        }
        for spans in chains.values_mut().flatten() {
            spans.sort_unstable();
            let mut merged = Vec::<ResidueSpan>::with_capacity(spans.len());
            for span in std::mem::take(spans) {
                if let Some(previous) = merged.last_mut()
                    && (upper_bound_matches(
                        span.start.number,
                        span.start.insertion.as_deref(),
                        &previous.end,
                    ) || previous.end.insertion.is_none()
                        && span.start.insertion.is_none()
                        && previous.end.number.checked_add(1) == Some(span.start.number))
                {
                    if upper_bound_precedes(&previous.end, &span.end) {
                        previous.end = span.end;
                    }
                } else {
                    merged.push(span);
                }
            }
            *spans = merged;
        }
        Ok(Self { chains })
    }

    pub(super) fn matches(&self, chain: &str, number: isize, insertion: Option<&str>) -> bool {
        if self.chains.is_empty() {
            return true;
        }
        match self.chains.get(chain) {
            Some(None) => true,
            Some(Some(spans)) => {
                let index = spans
                    .partition_point(|span| lower_bound_matches(number, insertion, &span.start));
                index > 0 && spans[index - 1].matches(number, insertion)
            }
            None => false,
        }
    }

    pub(super) fn validate_chains<'a>(
        &self,
        chains: impl Iterator<Item = &'a str>,
    ) -> ArpeggiaResult<()> {
        if self.chains.is_empty() {
            return Ok(());
        }
        let present = chains.collect::<BTreeSet<_>>();
        let unknown = self
            .chains
            .keys()
            .filter(|chain| !present.contains(chain.as_str()))
            .cloned()
            .collect::<Vec<_>>();
        if unknown.is_empty() {
            Ok(())
        } else {
            Err(ArpeggiaError::InvalidArgument(format!(
                "unknown chain identifiers: {}",
                unknown.join(",")
            )))
        }
    }
}

/// Validate both RMSD residue-selection grammars without structure I/O.
pub fn validate_rmsd_selections(
    superpose_residues: &str,
    rmsd_residues: &str,
) -> ArpeggiaResult<()> {
    ResidueSelector::parse(superpose_residues)
        .map_err(|error| selection_error("Superposition Selection", error))?;
    ResidueSelector::parse(rmsd_residues)
        .map_err(|error| selection_error("RMSD Selection", error))?;
    Ok(())
}

fn upper_bound_precedes(left: &ResidueBound, right: &ResidueBound) -> bool {
    left.number < right.number
        || left.number == right.number
            && match (&left.insertion, &right.insertion) {
                (None, _) => false,
                (Some(_), None) => true,
                (Some(left), Some(right)) => left < right,
            }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn residue_selector_supports_multichain_ranges() {
        let selector = ResidueSelector::parse("A:1-100,A:110-120,B,C:1,C:3,C:9-20").unwrap();
        assert!(selector.matches("A", 50, None));
        assert!(!selector.matches("A", 105, None));
        assert!(selector.matches("B", -999, None));
        assert!(selector.matches("C", 3, None));
        assert!(!selector.matches("C", 4, None));
        assert!(selector.matches("C", 12, None));
    }

    #[test]
    fn residue_selector_merges_overlapping_ranges() {
        let selector = ResidueSelector::parse("A:1-100,A:20-40,A:50-120,A:1-100").unwrap();
        assert_eq!(selector.chains["A"].as_ref().unwrap().len(), 1);
        assert!(selector.matches("A", 110, None));
        assert!(!selector.matches("A", 121, None));
    }

    #[test]
    fn residue_selector_merges_only_truly_adjacent_ranges() {
        let adjacent = ResidueSelector::parse("A:1-5,A:6-10").unwrap();
        assert_eq!(adjacent, ResidueSelector::parse("A:1-10").unwrap());

        let insertion_gap = ResidueSelector::parse("A:1-5A,A:6-10").unwrap();
        assert_eq!(insertion_gap.chains["A"].as_ref().unwrap().len(), 2);
    }

    #[test]
    fn residue_selector_supports_negative_numbers_and_insertions() {
        let selector = ResidueSelector::parse("A:-5--1,B:10A-20").unwrap();
        assert!(selector.matches("A", -3, None));
        assert!(!selector.matches("A", 1, None));
        assert!(!selector.matches("B", 10, None));
        assert!(selector.matches("B", 10, Some("A")));
        assert!(selector.matches("B", 20, Some("B")));
    }

    #[test]
    fn residue_selector_treats_bare_upper_bound_as_last_insertion() {
        let selector = ResidueSelector::parse("A:10A-10").unwrap();
        assert!(!selector.matches("A", 10, None));
        assert!(selector.matches("A", 10, Some("A")));
        assert!(selector.matches("A", 10, Some("B")));
        assert!(!selector.matches("A", 11, None));
    }

    #[test]
    fn residue_selector_rejects_malformed_and_unknown_chains() {
        for invalid in ["A:", "A:5-1", "A,,B", ":1"] {
            assert!(ResidueSelector::parse(invalid).is_err());
        }
        let selector = ResidueSelector::parse("missing").unwrap();
        assert!(selector.validate_chains(["A", "B"].into_iter()).is_err());
    }

    #[test]
    fn atom_subsets_include_caps_but_exclude_nonprotein_groups() {
        let input = std::env::temp_dir().join(format!(
            "arpeggia-rmsd-atom-subsets-{}.pdb",
            std::process::id()
        ));
        std::fs::write(
            &input,
            "HETATM    1  CH3 ACE A   0       0.000   0.000   0.000  1.00 20.00           C  \n\
             HETATM    2  H1  ACE A   0       0.000   1.000   0.000  1.00 20.00           H  \n\
             ATOM      3  N   ALA A   1       1.000   0.000   0.000  1.00 20.00           N  \n\
             ATOM      4  CA  ALA A   1       0.000   0.000   1.000  1.00 20.00           C  \n\
             HETATM    5  O   HOH A   2       3.000   3.000   3.000  1.00 20.00           O  \n\
             END\n",
        )
        .unwrap();
        let pdb = crate::load_model(input.to_str().unwrap()).unwrap().value;
        let selector = ResidueSelector::default();
        assert_eq!(
            select_coordinate_union(&pdb, 0, &selector, &selector, AtomSubset::Heavy)
                .unwrap()
                .keys
                .len(),
            3
        );
        assert_eq!(
            select_coordinate_union(&pdb, 0, &selector, &selector, AtomSubset::All)
                .unwrap()
                .keys
                .len(),
            4
        );
        assert_eq!(
            select_coordinate_union(&pdb, 0, &selector, &selector, AtomSubset::Backbone)
                .unwrap()
                .keys
                .len(),
            2
        );
    }

    #[test]
    fn coordinate_union_uses_three_contiguous_membership_blocks() {
        let input = std::env::temp_dir().join(format!(
            "arpeggia-rmsd-selection-union-{}.pdb",
            std::process::id()
        ));
        std::fs::write(
            &input,
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n\
             ATOM      2  CA  ALA B   1       1.000   0.000   0.000  1.00 20.00           C  \n\
             ATOM      3  CA  ALA C   1       0.000   1.000   0.000  1.00 20.00           C  \n\
             END\n",
        )
        .unwrap();
        let pdb = crate::load_model(input.to_str().unwrap()).unwrap().value;
        let superpose = ResidueSelector::parse("A,B").unwrap();
        let rmsd = ResidueSelector::parse("B,C").unwrap();
        let selected = select_coordinate_union(&pdb, 0, &superpose, &rmsd, AtomSubset::Ca).unwrap();
        assert_eq!(selected.superpose_end, 2);
        assert_eq!(selected.rmsd_start, 1);
        assert_eq!(
            selected.coordinates.len() * std::mem::size_of::<[f64; 3]>(),
            72
        );
        assert_eq!(
            selected
                .keys
                .iter()
                .map(|identity| identity.chain.as_str())
                .collect::<Vec<_>>(),
            ["A", "B", "C"]
        );

        let missing = ResidueSelector::parse("missing").unwrap();
        assert!(matches!(
            select_coordinate_union(&pdb, 0, &missing, &rmsd, AtomSubset::Ca),
            Err(ArpeggiaError::InvalidArgument(message))
                if message.starts_with("Superposition Selection:")
        ));
        assert!(matches!(
            select_coordinate_union(&pdb, 0, &superpose, &missing, AtomSubset::Ca),
            Err(ArpeggiaError::InvalidArgument(message))
                if message.starts_with("RMSD Selection:")
        ));
    }

    #[test]
    fn heavy_subset_infers_digit_leading_hydrogen_names() {
        let hydrogen = pdbtbx::Atom::new(false, 1, "1HB", 0.0, 0.0, 0.0, 1.0, 20.0, "", 0).unwrap();
        assert!(hydrogen.element().is_none());
        assert!(!atom_in_subset(
            hydrogen.name(),
            hydrogen.element(),
            AtomSubset::Heavy
        ));
        assert!(!atom_in_subset("D1", None, AtomSubset::Heavy));
        assert!(atom_in_subset("CB", None, AtomSubset::Heavy));
    }

    #[test]
    fn exact_correspondence_rejects_missing_or_different_atoms() {
        let atom = AtomIdentity {
            chain: "A".into(),
            residue_number: 1,
            insertion: String::new(),
            residue_name: "ALA".into(),
            atom_name: "CA".into(),
        };
        let mut different = atom.clone();
        different.atom_name = "N".into();
        assert!(validate_correspondence(std::slice::from_ref(&atom), &[]).is_err());
        assert!(validate_correspondence(&[atom], &[different]).is_err());
    }
}
