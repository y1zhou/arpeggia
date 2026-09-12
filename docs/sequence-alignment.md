# Pairwise sequence alignment

## Align two sequences

```python
import arpeggia

alignment = arpeggia.align_seqs("GGACDEFGHIKGG", "ACDEFGHIK", mode="semi-global")
print(alignment.score, alignment.identity_alignment, alignment.edit_distance)
print(alignment.reference_span, alignment.query_span)  # (2, 11), (0, 9)
print(alignment.aligned_reference, alignment.aligned_query)
print(alignment.operations)  # space match, + insertion, - deletion, : similar, x other
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
places and require opening ≥ extension ≥ 0.01; smaller costs are rejected, never
rounded to zero. Matrix and costs use the same integer scale; unrepresentable
scores fail explicitly.

Lowercase is normalized; standard amino acids and `B/Z/X/U/O` are accepted.
`U/O` score as `C/K` with a warning but retain their original symbols for identity
and edit distance. Gaps, stop symbols, whitespace, unsupported characters, and
empty inputs are invalid. One deterministic optimum is returned; equal scores
do not establish a unique biological correspondence.

`SeqAlignment` has read-only Python attributes. Spans are zero-based and
half-open. `aligned_reference` and `aligned_query` contain the scored alignment
with `-` for gaps; `operations` describes reference-to-query changes. These three
ASCII strings have equal lengths. Clipped terminal segments remain in the original
inputs, outside the spans and aligned strings. Input indices can be recovered
by counting nongap residues from each span start.

| Attribute | Meaning |
| --- | --- |
| `score` | Optimal BLOSUM62/affine score, in original matrix units |
| `matches`, `mismatches`, `paired_residues` | Identical, substituted, and total nongap pairs |
| `alignment_length`, `shorter_length` | Column count including gaps; shorter full input length |
| `identity_alignment`, `identity_shorter` | Matches divided by each named denominator |
| `coverage_alignment`, `coverage_shorter` | Nongap pairs divided by each named denominator |
| `gap_residues`, `gap_runs` | Residues opposite gaps; contiguous gap runs |
| `edit_distance` | Minimum full-input single-residue substitutions, insertions, and deletions |

A local alignment without a positive match has score zero, empty aligned strings,
`None` alignment-length ratios, and zero shorter-input ratios. Edit distance
still compares the complete inputs. JSON represents undefined ratios as `null`.

## Display an alignment

Python `SeqAlignment` displays statistics, reference and query rows, and an
unlabeled operation row. The CLI uses the same layout, including per-chain
alignments from `rmsd --align-seqs`.

```python
alignment = arpeggia.align_seqs(
    "ACDEFGHIKLMN", "ACDFGHIKLMN",
    reference_name="Wild type", query_name="Mutant",
)
print(alignment)  # terminal width and automatic color
print(alignment.format(width=60, color="never", rulers=False))
```

```bash
arpeggia align-seqs ACDEFGHIKLMN ACDFGHIKLMN \
  --reference-name "Wild type" --query-name "Mutant" --width 60 --no-rulers
```

`reference_name` and `query_name` are optional keyword-only display names,
defaulting to "Reference" and "Query". Both are read-only result fields and appear
in JSON. They do not change sequence data, scores, or correspondence. Labels
are padded to their visible widths; control characters are escaped for display.
RMSD chain alignments use names such as "Reference A" and "Query H".

Operations describe reference-to-query changes: deletions are red (`-`),
insertions green (`+`), similar substitutions blue (`:`), other substitutions
yellow (`x`), and matches uncolored (space).
Both sequence cells and the operation marker share the color. Unaligned tails
are gray with blank operation cells; they do not contribute to alignment
statistics. Terminal gaps in global alignment remain scored edits.

Similar substitutions are nonidentical pairs with a strictly positive BLOSUM62
score, including ambiguous `B/Z` residues and `U/O` scoring aliases. Zero scores remain `x`.
Identity takes precedence: identical `X/X` is a match despite its negative score.
Both substitution categories count toward `mismatches`; the distinction changes
neither identity nor edit distance and does not establish chemical equivalence
for structural atom pairing. Stored operations and JSON use the same ASCII markers.

Each block shows one-based input start/end positions and separate rulers above
each sequence, with every tenth position right-aligned to its residue. Gaps and
padding do not advance positions. These are sequence positions, including for
structural alignments; author residue numbers are in `residue_pairs`. API spans
remain zero-based and half-open.

Width includes labels and position numbers. Output wraps without truncation,
using terminal width or 80 columns if unavailable; widths too small for labels
and one residue fail. `--no-rulers` and `rulers=False` hide tick rows while
retaining endpoints. `--color` and `color` accept `auto` (default), `always`, or
`never`. Auto respects terminal capability and `NO_COLOR`; explicit `always`
forces ANSI color. Fields and JSON always contain plain data.

## Establish structural correspondence

Residue selectors use comma-separated chains and inclusive author-number ranges,
e.g. `A:1-100,A:110-120,B`. Repeat the chain in each clause. Empty selectors
independently mean all eligible residues; an empty evaluation selection never
inherits the fitting selection. See the [selection guide](https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md#rmsd)
for insertion codes, negative numbers, and atom populations.

```python
result = arpeggia.rmsd(
    "reference.cif", "query.cif",
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
arpeggia rmsd reference.cif query.cif \
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
selection and assign unique query partners. Without a map, Arpeggia scores
all relevant reference/query chain pairs semi-globally, consuming shorter against
longer, and maximizes the total score of a one-to-one assignment. All inferred
pair scores must be positive. An absent complete assignment or tied optimum
requires an explicit map; unused query chains are allowed. Positive scores
alone do not establish homology: inspect the returned identity and coverage.

Final residue alignments follow `alignment_mode` (default global), independently
of the semi-global scoring used to infer chain partners. Reference remains first,
query second, so final semi-global alignment consumes the full query chain.

## Refinement and result migration

`refine_cycles=0` performs one initial fit without rejection. Each additional
cycle rejects fitting atom pairs farther than `refine_cutoff × current fitting
RMSD` (default multiplier 2), then refits survivors. Stop when unchanged or
numerically exact; fail if survivors cannot define a non-collinear fit.
Rejection is permanent and atom-wise. It also works without sequence alignment,
using the existing exact correspondence; it does not infer structural matches
as PyMOL `super` does.

`rmsd()` returns a read-only `RmsdResult`. When migrating from the scalar API, use
`result.rmsd`; Rust callers use `get_rmsd(reference, query, &RmsdOptions)` and
read `analysis.value.rmsd`. Python callers using the old `mobile=` keyword
should use `query=` for the second structure.

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

`pairwise_rmsd` and clustering require exact correspondence; sequence-aware
collection comparisons remain deferred. [Antibody numbering](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md)
provides a separate numbered-position correspondence. See
[ADR 0009](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0009-separate-sequence-correspondence-from-rmsd-evaluation.md) for
rationale and [validation](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/sequence-alignment.md) for reference checks.
