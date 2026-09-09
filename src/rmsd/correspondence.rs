//! Observed-chain alignment and coordinate correspondence before atom selection.
use super::*;
use pdbtbx::{Atom, Chain, Residue};

/// Corresponding coordinate residues, preserving author numbering and names.
#[derive(Clone, Debug, Serialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object, module = "arpeggia")
)]
pub struct ResiduePair {
    /// Reference author residue number.
    pub reference_number: isize,
    /// Reference insertion code (empty when absent).
    pub reference_insertion: String,
    /// Reference chemical residue name.
    pub reference_name: String,
    /// Query author residue number.
    pub query_number: isize,
    /// Query insertion code.
    pub query_insertion: String,
    /// Query chemical residue name.
    pub query_name: String,
}

/// One paired chain and its observed-sequence/residue correspondence.
#[derive(Clone, Debug, Serialize)]
#[cfg_attr(
    feature = "python",
    pyo3::pyclass(frozen, get_all, skip_from_py_object, module = "arpeggia")
)]
pub struct ChainAlignment {
    /// Reference chain ID.
    pub reference_chain: String,
    /// Query chain ID.
    pub query_chain: String,
    /// Final sequence alignment; None for exact atom correspondence.
    pub alignment: Option<crate::SeqAlignment>,
    /// Nongap coordinate-residue pairs in the final sequence alignment.
    pub residue_pairs: Vec<ResiduePair>,
}

pub(super) fn exact_chains(keys: &[AtomIdentity]) -> Vec<ChainAlignment> {
    let residues: BTreeSet<_> = keys
        .iter()
        .map(|k| (&k.chain, k.residue_number, &k.insertion, &k.residue_name))
        .collect();
    let mut chains: BTreeMap<String, ChainAlignment> = BTreeMap::new();
    for (chain, number, insertion, name) in residues {
        chains
            .entry(chain.clone())
            .or_insert_with(|| ChainAlignment {
                reference_chain: chain.clone(),
                query_chain: chain.clone(),
                alignment: None,
                residue_pairs: Vec::new(),
            })
            .residue_pairs
            .push(ResiduePair {
                reference_number: number,
                reference_insertion: insertion.clone(),
                reference_name: name.clone(),
                query_number: number,
                query_insertion: insertion.clone(),
                query_name: name.clone(),
            });
    }
    chains.into_values().collect()
}

struct ObservedChain<'a> {
    chain: &'a Chain,
    residues: Vec<&'a Residue>,
    sequence: String,
    encoded: Vec<u8>,
}
fn observed_chains(pdb: &PDB, model: usize) -> ArpeggiaResult<Vec<ObservedChain<'_>>> {
    let mut ids = BTreeSet::new();
    let mut result = Vec::new();
    for chain in selected_model(pdb, model)?.chains() {
        if !ids.insert(chain.id()) {
            return Err(ArpeggiaError::Calculation(format!(
                "ambiguous duplicate chain ID {:?}",
                chain.id()
            )));
        }
        let mut residues = Vec::new();
        let mut sequence = String::new();
        let mut seen = BTreeSet::new();
        for residue in chain.residues() {
            if let Some(letter) = residue.name().and_then(one_letter_code) {
                if !seen.insert((residue.serial_number(), residue.insertion_code())) {
                    return Err(ArpeggiaError::Calculation(format!(
                        "duplicate residue identity in chain {}",
                        chain.id()
                    )));
                }
                residues.push(residue);
                sequence.push_str(letter);
            }
        }
        if !residues.is_empty() {
            let encoded = crate::seq_alignment::encode(&sequence)?;
            result.push(ObservedChain {
                chain,
                residues,
                sequence,
                encoded,
            });
        }
    }
    result.sort_by(|a, b| a.chain.id().cmp(b.chain.id()));
    Ok(result)
}

// Rectangular Hungarian assignment in O(n^2 m), using integer scores. No current
// dependency supplies assignment; this avoids importing a second matrix stack.
fn assignment(scores: &[Vec<i32>], forbidden: Option<(usize, usize)>) -> Option<(i64, Vec<usize>)> {
    let n = scores.len();
    let m = scores.first()?.len();
    if n > m {
        return None;
    }
    const INF: i64 = i64::MAX / 8;
    let mut u = vec![0i64; n + 1];
    let mut v = vec![0i64; m + 1];
    let mut p = vec![0usize; m + 1];
    let mut way = vec![0usize; m + 1];
    for i in 1..=n {
        p[0] = i;
        let mut j0 = 0;
        let mut minv = vec![INF; m + 1];
        let mut used = vec![false; m + 1];
        loop {
            used[j0] = true;
            let i0 = p[j0];
            let mut delta = INF;
            let mut j1 = 0;
            for j in 1..=m {
                if used[j] {
                    continue;
                }
                let score = scores[i0 - 1][j - 1];
                if score > 0 && forbidden != Some((i0 - 1, j - 1)) {
                    let cost = -i64::from(score) - u[i0] - v[j];
                    if cost < minv[j] {
                        minv[j] = cost;
                        way[j] = j0;
                    }
                }
                if minv[j] < delta {
                    delta = minv[j];
                    j1 = j;
                }
            }
            if delta == INF {
                return None;
            }
            for j in 0..=m {
                if used[j] {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else if minv[j] != INF {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
            if p[j0] == 0 {
                break;
            }
        }
        loop {
            let j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
            if j0 == 0 {
                break;
            }
        }
    }
    let mut map = vec![0; n];
    for j in 1..=m {
        if p[j] > 0 {
            map[p[j] - 1] = j - 1;
        }
    }
    let score = map
        .iter()
        .enumerate()
        .map(|(i, j)| i64::from(scores[i][*j]))
        .sum();
    Some((score, map))
}

fn chain_mapping(
    first: &[ObservedChain<'_>],
    second: &[ObservedChain<'_>],
    options: &RmsdOptions,
) -> ArpeggiaResult<Vec<usize>> {
    if !options.chain_map.is_empty() {
        if options.chain_map.len() != first.len() {
            return Err(ArpeggiaError::InvalidArgument(
                "chain_map must cover exactly the selected reference chains".into(),
            ));
        }
        let mut used = BTreeSet::new();
        return first
            .iter()
            .map(|chain| {
                let id = options.chain_map.get(chain.chain.id()).ok_or_else(|| {
                    ArpeggiaError::InvalidArgument(format!(
                        "chain_map is missing reference chain {}",
                        chain.chain.id()
                    ))
                })?;
                if !used.insert(id) {
                    return Err(ArpeggiaError::InvalidArgument(
                        "chain_map query partners must be unique".into(),
                    ));
                }
                second
                    .iter()
                    .position(|chain| chain.chain.id() == id)
                    .ok_or_else(|| {
                        ArpeggiaError::InvalidArgument(format!(
                            "unknown or nonprotein query chain {id}"
                        ))
                    })
            })
            .collect();
    }
    let scoring = crate::seq_alignment::scoring(&options.alignment)?;
    let mut scratch = hyalite::PairScratch::new();
    let mut scores = Vec::new();
    for a in first {
        let mut row = Vec::new();
        for b in second {
            let (long, short) = if a.encoded.len() >= b.encoded.len() {
                (&a.encoded, &b.encoded)
            } else {
                (&b.encoded, &a.encoded)
            };
            row.push(
                hyalite::align_pair_with(
                    &mut scratch,
                    long,
                    short,
                    &scoring,
                    hyalite::Mode::Shw,
                    hyalite::SearchType::Score,
                )
                .map_err(crate::seq_alignment::backend_error)?
                .score,
            );
        }
        scores.push(row);
    }
    let (best, map) = assignment(&scores, None).ok_or_else(|| {
        ArpeggiaError::Calculation(
            "no complete positive-score chain assignment; supply an explicit chain_map".into(),
        )
    })?;
    // Any different optimum excludes at least one edge of this optimum. Re-solving
    // with each selected edge forbidden detects ties without enumerating maps.
    for (i, j) in map.iter().enumerate() {
        if assignment(&scores, Some((i, *j))).is_some_and(|(score, _)| score == best) {
            return Err(ArpeggiaError::Calculation(
                "ambiguous optimal chain assignment; supply an explicit chain_map".into(),
            ));
        }
    }
    Ok(map)
}

fn residue_pair(a: &Residue, b: &Residue) -> ResiduePair {
    ResiduePair {
        reference_number: a.serial_number(),
        reference_insertion: a.insertion_code().unwrap_or("").into(),
        reference_name: a.name().unwrap_or("").into(),
        query_number: b.serial_number(),
        query_insertion: b.insertion_code().unwrap_or("").into(),
        query_name: b.name().unwrap_or("").into(),
    }
}

fn selected_atoms(residue: &Residue, subset: AtomSubset) -> ArpeggiaResult<BTreeMap<&str, &Atom>> {
    let mut result = BTreeMap::new();
    for atom in residue
        .atoms()
        .filter(|a| atom_in_subset(a.name(), a.element(), subset))
    {
        if result.insert(atom.name(), atom).is_some() {
            return Err(ArpeggiaError::Calculation(format!(
                "duplicate selected atom {} in residue {}",
                atom.name(),
                residue.serial_number()
            )));
        }
    }
    Ok(result)
}

pub(super) fn aligned_coordinates(
    reference: &PDB,
    query: &PDB,
    options: &RmsdOptions,
    fit: &ResidueSelector,
    eval: &ResidueSelector,
    warnings: &mut Vec<AnalysisWarning>,
) -> ArpeggiaResult<(
    SelectedCoordinateUnion,
    SelectedCoordinateUnion,
    Vec<ChainAlignment>,
)> {
    let model = selected_model(reference, options.model_num)?;
    fit.validate_chains(model.chains().map(|c| c.id()))?;
    eval.validate_chains(model.chains().map(|c| c.id()))?;
    let mut first = observed_chains(reference, options.model_num)?;
    first.retain(|c| {
        c.residues.iter().any(|r| {
            fit.matches(c.chain.id(), r.serial_number(), r.insertion_code())
                || eval.matches(c.chain.id(), r.serial_number(), r.insertion_code())
        })
    });
    if first.is_empty() {
        return Err(ArpeggiaError::Calculation(
            "no observed protein residues in selections".into(),
        ));
    }
    let second = observed_chains(query, options.model_num)?;
    let map = chain_mapping(&first, &second, options)?;
    let mut selected = Vec::new();
    let mut chains = Vec::new();
    let mut omitted = 0usize;
    for (a, j) in first.iter().zip(map) {
        let b = &second[j];
        let analysis = crate::align_seqs(&a.sequence, &b.sequence, &options.alignment)?;
        warnings.extend(analysis.warnings);
        let alignment = crate::SeqAlignment {
            reference_name: format!("Reference {}", a.chain.id()),
            query_name: format!("Query {}", b.chain.id()),
            ..analysis.value
        };
        let mut residue_pairs = Vec::new();
        for (i, j) in alignment.columns() {
            let (Some(i), Some(j)) = (i, j) else {
                continue;
            };
            let left = a.residues[i];
            let right = b.residues[j];
            residue_pairs.push(residue_pair(left, right));
            let membership = match (
                fit.matches(a.chain.id(), left.serial_number(), left.insertion_code()),
                eval.matches(a.chain.id(), left.serial_number(), left.insertion_code()),
            ) {
                (true, true) => SelectionMembership::Both,
                (true, false) => SelectionMembership::SuperposeOnly,
                (false, true) => SelectionMembership::RmsdOnly,
                (false, false) => continue,
            };
            let left_atoms = selected_atoms(left, options.atoms)?;
            let right_atoms = selected_atoms(right, options.atoms)?;
            let equivalent = a.sequence.as_bytes()[i] == b.sequence.as_bytes()[j];
            let mut paired = 0;
            for (name, atom) in &left_atoms {
                let Some(other) = right_atoms.get(name) else {
                    continue;
                };
                if atom.element() != other.element()
                    || (!equivalent && !atom_in_subset(name, atom.element(), AtomSubset::Backbone))
                {
                    continue;
                }
                let coords = [atom.x(), atom.y(), atom.z()];
                let other_coords = [other.x(), other.y(), other.z()];
                if coords.iter().chain(&other_coords).any(|x| !x.is_finite()) {
                    return Err(ArpeggiaError::Calculation(
                        "nonfinite paired atom coordinates".into(),
                    ));
                }
                selected.push((
                    membership,
                    AtomIdentity {
                        chain: a.chain.id().into(),
                        residue_number: left.serial_number(),
                        insertion: left.insertion_code().unwrap_or("").into(),
                        residue_name: left.name().unwrap_or("").into(),
                        atom_name: (*name).into(),
                    },
                    coords,
                    other_coords,
                ));
                paired += 1;
            }
            omitted += left_atoms.len() + right_atoms.len() - 2 * paired;
        }
        // Include selected reference residues excluded by gaps or terminal clipping.
        let paired_indices: BTreeSet<_> =
            alignment.columns().filter_map(|(i, j)| j.and(i)).collect();
        for (i, r) in a.residues.iter().enumerate() {
            if !paired_indices.contains(&i)
                && (fit.matches(a.chain.id(), r.serial_number(), r.insertion_code())
                    || eval.matches(a.chain.id(), r.serial_number(), r.insertion_code()))
            {
                omitted += selected_atoms(r, options.atoms)?.len();
            }
        }
        chains.push(ChainAlignment {
            reference_chain: a.chain.id().into(),
            query_chain: b.chain.id().into(),
            alignment: Some(alignment),
            residue_pairs,
        });
    }
    if omitted > 0 {
        warnings.push(AnalysisWarning::new(WarningCode::IncompleteCorrespondence,format!("sequence correspondence omitted {omitted} selected atom endpoints; inspect alignment coverage")));
    }
    for pdb in [reference, query] {
        if options.model_num == 0 && pdb.model_count() > 1 {
            warnings.push(AnalysisWarning::new(
                WarningCode::ModelSelected,
                format!(
                    "selected first model {}",
                    selected_model(pdb, 0)?.serial_number()
                ),
            ));
        }
    }
    selected.sort_by(|a, b| (&a.0, &a.1).cmp(&(&b.0, &b.1)));
    let fit_end = selected.partition_point(|row| row.0 != SelectionMembership::RmsdOnly);
    let eval_start = selected.partition_point(|row| row.0 == SelectionMembership::SuperposeOnly);
    let keys: Vec<_> = selected.iter().map(|row| row.1.clone()).collect();
    let first = SelectedCoordinateUnion {
        keys: keys.clone(),
        coordinates: selected.iter().map(|row| row.2).collect(),
        superpose_end: fit_end,
        rmsd_start: eval_start,
        warnings: Vec::new(),
    };
    let second = SelectedCoordinateUnion {
        keys,
        coordinates: selected.iter().map(|row| row.3).collect(),
        superpose_end: fit_end,
        rmsd_start: eval_start,
        warnings: Vec::new(),
    };
    Ok((first, second, chains))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn structure(chains: &[(&str, &[&str], isize)]) -> PDB {
        let mut model = pdbtbx::Model::new(1);
        let mut serial = 1;
        for (c, names, offset) in chains {
            for (i, name) in names.iter().enumerate() {
                let atom = Atom::new(
                    false,
                    serial,
                    "CA",
                    (i % 2) as f64,
                    (i / 2) as f64,
                    0.0,
                    1.0,
                    20.0,
                    "C",
                    0,
                )
                .unwrap();
                model.add_atom(atom, c, (*offset + i as isize, None), (*name, None));
                serial += 1;
            }
        }
        let mut pdb = PDB::new();
        pdb.add_model(model);
        pdb
    }
    #[test]
    fn optimal_assignment_is_not_greedy_and_rejects_missing_edges() {
        let scores = vec![vec![100, 99], vec![98, 1]];
        assert_eq!(assignment(&scores, None), Some((197, vec![1, 0])));
        assert_eq!(assignment(&[vec![0, 3], vec![-1, 2]], None), None);
        assert_eq!(assignment(&[vec![3, 2, 1]], None), Some((3, vec![0])));
        assert_eq!(assignment(&[vec![3], vec![2]], None), None);
    }
    #[test]
    fn assignment_matches_exhaustive_small_matrices() {
        // Exhaustive score alphabets include prohibited, tied, and negative edges.
        for bits in 0usize..4096 {
            let mut x = bits;
            let mut scores = vec![vec![0; 3]; 2];
            for row in &mut scores {
                for v in row {
                    *v = (x % 4) as i32 - 1;
                    x /= 4;
                }
            }
            let mut best = None;
            for i in 0..3 {
                for j in 0..3 {
                    if i != j && scores[0][i] > 0 && scores[1][j] > 0 {
                        let score = i64::from(scores[0][i] + scores[1][j]);
                        best = Some(best.map_or(score, |v: i64| v.max(score)));
                    }
                }
            }
            assert_eq!(assignment(&scores, None).map(|(v, _)| v), best);
        }
    }
    #[test]
    fn renumbered_mutant_uses_reference_selections() {
        let reference = structure(&[("A", &["ALA", "CYS", "ASP", "GLU"], 1)]);
        let query = structure(&[("H", &["ALA", "CYS", "ASN", "GLU"], 101)]);
        assert!(get_rmsd(reference.clone(), query.clone(), &RmsdOptions::default()).is_err());
        let result = get_rmsd(
            reference,
            query,
            &RmsdOptions {
                align_seqs: true,
                superpose_residues: "A:1-3".into(),
                rmsd_residues: "A:4".into(),
                ..Default::default()
            },
        )
        .unwrap()
        .value;
        assert!(result.rmsd < 1e-12);
        assert_eq!((result.initial_fit_atoms, result.evaluation_atoms), (3, 1));
        let chain = &result.chain_alignments[0];
        assert_eq!(chain.query_chain, "H");
        assert_eq!(chain.residue_pairs[2].query_number, 103);
        assert_eq!(chain.alignment.as_ref().unwrap().mismatches, 1);
    }
    #[test]
    fn inferred_chains_are_unique_and_explicit_maps_override_ambiguity() {
        let first = structure(&[
            ("A", &["ALA", "CYS", "ASP"], 1),
            ("B", &["TRP", "TYR", "PHE"], 1),
        ]);
        let second = structure(&[
            ("X", &["TRP", "TYR", "PHE"], 100),
            ("Y", &["ALA", "CYS", "ASP"], 200),
            ("Z", &["GLY", "GLY", "GLY"], 1),
        ]);
        let result = get_rmsd(
            first,
            second,
            &RmsdOptions {
                align_seqs: true,
                ..Default::default()
            },
        )
        .unwrap()
        .value;
        assert_eq!(
            result
                .chain_alignments
                .iter()
                .map(|c| (c.reference_chain.as_str(), c.query_chain.as_str()))
                .collect::<Vec<_>>(),
            vec![("A", "Y"), ("B", "X")]
        );
        let first = structure(&[("A", &["ALA", "CYS", "ASP"], 1)]);
        let second = structure(&[
            ("X", &["ALA", "CYS", "ASP"], 100),
            ("Y", &["ALA", "CYS", "ASP"], 200),
        ]);
        let mut opts = RmsdOptions {
            align_seqs: true,
            ..Default::default()
        };
        assert!(
            get_rmsd(first.clone(), second.clone(), &opts)
                .unwrap_err()
                .to_string()
                .contains("ambiguous")
        );
        opts.chain_map.insert("A".into(), "Y".into());
        assert_eq!(
            get_rmsd(first, second, &opts)
                .unwrap()
                .value
                .chain_alignments[0]
                .query_chain,
            "Y"
        );
    }
    #[test]
    fn gaps_and_substitution_sidechains_are_omitted_with_diagnostics() {
        let first = structure(&[("A", &["ALA", "CYS", "ASP", "GLU", "PHE"], 1)]);
        let second = structure(&[("H", &["ALA", "CYS", "GLU", "PHE"], 101)]);
        let analysis = get_rmsd(
            first,
            second,
            &RmsdOptions {
                align_seqs: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(analysis.value.evaluation_atoms, 4);
        assert!(
            analysis
                .warnings
                .iter()
                .any(|w| w.code == WarningCode::IncompleteCorrespondence)
        );
        let mut first = structure(&[("A", &["ALA", "ASP", "GLY"], 1)]);
        let mut second = structure(&[("H", &["ALA", "LEU", "GLY"], 101)]);
        for (pdb, c, number, name) in [(&mut first, "A", 2, "ASP"), (&mut second, "H", 102, "LEU")]
        {
            pdb.models_mut().next().unwrap().add_atom(
                Atom::new(false, 10, "CG", 1.0, 0.0, 1.0, 1.0, 20.0, "C", 0).unwrap(),
                c,
                (number, None),
                (name, None),
            );
        }
        let analysis = get_rmsd(
            first,
            second,
            &RmsdOptions {
                align_seqs: true,
                atoms: AtomSubset::Heavy,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(analysis.value.evaluation_atoms, 3);
        assert!(
            analysis
                .warnings
                .iter()
                .any(|w| w.code == WarningCode::IncompleteCorrespondence)
        );
    }
    #[test]
    fn refinement_preserves_all_evaluation_pairs_and_can_fail() {
        let directions = [
            [1., 0., 0.],
            [0., 1., 0.],
            [0., 0., 1.],
            [1., 1., 0.],
            [1., 0., 1.],
            [0., 1., 1.],
            [1., -1., 0.],
            [1., 0., -1.],
            [0., 1., -1.],
            [1., 1., 1.],
        ];
        let mut reference = pdbtbx::Model::new(1);
        let mut query = pdbtbx::Model::new(1);
        let mut serial = 1;
        for (i, direction) in directions.iter().enumerate() {
            let norm: f64 = direction.iter().map(|v| v * v).sum::<f64>().sqrt();
            let delta = if i < 8 {
                0.1
            } else if i == 8 {
                1.0
            } else {
                6.0
            };
            for sign in [-1., 1.] {
                for (model, radius) in [(&mut reference, 10.0), (&mut query, 10.0 + delta)] {
                    let p = direction.map(|v| v / norm * radius * sign);
                    model.add_atom(
                        Atom::new(false, serial, "CA", p[0], p[1], p[2], 1.0, 20.0, "C", 0)
                            .unwrap(),
                        "A",
                        (serial as isize, None),
                        ("ALA", None),
                    );
                }
                serial += 1;
            }
        }
        let mut a = PDB::new();
        a.add_model(reference);
        let mut b = PDB::new();
        b.add_model(query);
        let result = get_rmsd(
            a.clone(),
            b.clone(),
            &RmsdOptions {
                refine_cycles: 10,
                ..Default::default()
            },
        )
        .unwrap()
        .value;
        assert_eq!(result.evaluation_atoms, 20);
        assert_eq!(result.retained_fit_atoms, 16);
        assert_eq!(result.cycles, 3);
        assert!((result.core_rmsd - 0.1).abs() < 1e-10);
        assert!((result.rmsd - (74.16_f64 / 20.).sqrt()).abs() < 1e-10);
        let fail = get_rmsd(
            a,
            b,
            &RmsdOptions {
                refine_cycles: 1,
                refine_cutoff: 0.001,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(fail.to_string().contains("refinement retained 0"));
    }
}
