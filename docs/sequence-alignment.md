# Pairwise sequence alignment and sequence-aware RMSD

## Align two sequences

```python
import arpeggia

alignment = arpeggia.align_seqs("GGACDEFGHIKGG", "ACDEFGHIK", mode="semi-global")
print(alignment.score, alignment.identity_alignment, alignment.edit_distance)
print(alignment.reference_span, alignment.mobile_span)  # (2, 11), (0, 9)
print(alignment.columns[0])  # (2, 0): zero-based input indices
```

```bash
arpeggia align-seqs GGACDEFGHIKGG ACDEFGHIK --mode semi-global --json
```

Inputs are two unaligned amino-acid strings. `global` aligns both end to end;
`local` finds the best-scoring subsequences; `semi-global` consumes the complete
**second** sequence with free terminal overhangs of the first. Global is the
default. The CLI accepts `--alignment-mode` and its alias `--mode`.

BLOSUM62 scoring uses positive `gap_open=10` and `gap_extend=0.5`, charging
`gap_open + (length - 1) * gap_extend` per gap. Costs accept at most two decimal
places and require opening ≥ extension. Matrix and costs are scaled together
for exact integer optimization. Unrepresentable scores fail explicitly.

Lowercase is normalized; standard amino acids and `B/Z/X/U/O` are accepted.
`U/O` score as `C/K` with a warning but retain their original symbols for identity
and edit distance. Gaps, stop symbols, whitespace, unsupported characters, and
empty inputs are invalid. One deterministic optimum is returned; equal scores
do not establish a unique biological correspondence.

`SeqAlignment` has read-only Python attributes. Spans are zero-based and
half-open; `columns` contains an index pair per alignment column, with `None`
for gaps. Clipped terminal segments are outside the spans and columns.

| Attribute | Meaning |
| --- | --- |
| `score` | Optimal BLOSUM62/affine score, in original matrix units |
| `matches`, `mismatches`, `paired_residues` | Identical, substituted, and total nongap pairs |
| `alignment_length`, `shorter_length` | Column count including gaps; shorter full input length |
| `identity_alignment`, `identity_shorter` | Matches divided by each named denominator |
| `coverage_alignment`, `coverage_shorter` | Nongap pairs divided by each named denominator |
| `gap_residues`, `gap_runs` | Residues opposite gaps; contiguous gap runs |
| `edit_distance` | Minimum full-input single-residue substitutions, insertions, and deletions |

A local alignment without a positive match has score zero, empty columns,
`None` alignment-length ratios, and zero shorter-input ratios. Edit distance
still compares the complete inputs. JSON represents undefined ratios as `null`.

## Establish structural correspondence

```python
result = arpeggia.rmsd(
    "reference.cif", "mobile.cif",
    align_seqs=True,
    chain_map={"A": "H", "B": "L"},
    superpose_residues="A:1-100",
    rmsd_residues="B",
    atoms="ca",
    refine_cycles=2,
)
print(result.rmsd, result.core_rmsd, result.retained_fit_atoms)
```

```bash
arpeggia rmsd reference.cif mobile.cif \
  --align-seqs --chain-map A=H --chain-map B=L \
  --superpose-residues A:1-100 --rmsd-residues B \
  --refine-cycles 2 --json
```

With alignment enabled, both residue selectors address reference author
numbering. Complete observed chain sequences are aligned before applying residue
ranges and atom selection. Declared sequences and residues without coordinates
cannot supply atom pairs. Nongap substitutions remain corresponding residues.
Backbone atoms can pair across substitutions; side-chain pairing requires the
same normalized amino-acid type, atom name, and element. Missing counterparts
are omitted with diagnostics; sequence alignment does not reconstruct atoms.

An explicit chain map must cover exactly the reference chains used by either
selection and assign unique mobile partners. Without a map, Arpeggia scores
all relevant reference/mobile chain pairs semi-globally, consuming shorter against
longer, and maximizes the total score of a one-to-one assignment. All inferred
pair scores must be positive. An absent complete assignment or tied optimum
requires an explicit map; unused mobile chains are allowed. Positive scores
alone do not establish homology: inspect the returned identity and coverage.

Final residue alignments follow `alignment_mode` (default global), independently
of the semi-global scoring used to infer chain partners. Reference remains first,
mobile second, so final semi-global alignment consumes the full mobile chain.

## Refinement and result migration

`refine_cycles=0` performs one initial fit without rejection. Each additional
cycle rejects fitting atom pairs farther than `refine_cutoff × current fitting
RMSD` (default multiplier 2), then refits survivors. Stop when unchanged or
numerically exact; fail if survivors cannot define a non-collinear fit.
Rejection is permanent and atom-wise. It also works without sequence alignment,
using the existing exact correspondence; it does not infer structural matches
as PyMOL `super` does.

`rmsd()` now returns a read-only `RmsdResult`. Replace scalar uses with
`result.rmsd`; Rust callers use `get_rmsd(reference, mobile, &RmsdOptions)` and
read `analysis.value.rmsd`.

- `rmsd`: all mapped evaluation pairs under the final retained-pair transform.
- `core_rmsd`: retained fitting pairs under that transform.
- `initial_rmsd`: fitting RMSD before rejection.
- `initial_fit_atoms`, `retained_fit_atoms`, `evaluation_atoms`: atom-pair counts.
- `cycles`: rejection inspections performed, including a final unchanged pass.
- `chain_alignments`: paired chain IDs, optional sequence results, and coordinate
  residue pairs preserving author numbers, insertion codes, and residue names.

Pairs rejected from fitting remain in evaluation if selected. The evaluation
region is never refitted. Results include parameters and selected model serials,
without per-atom residuals or transformations. The CLI prints a detailed summary;
`--json` provides the structured result.

`pairwise_rmsd` and clustering still require exact correspondence. Sequence-aware
collection comparisons and antibody numbering are follow-ups. See
[ADR 0009](adr/0009-separate-sequence-correspondence-from-rmsd-evaluation.md) for
rationale and [validation](benchmarks/sequence-alignment.md) for reference checks.
