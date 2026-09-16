# Antibody numbering

`number_antibody()` identifies one antibody variable domain, numbers its
residues, assigns framework/CDR regions, and finds the closest bundled V and J
references. `align_antibodies()` compares numbered positions across antibodies.
Python aligns `NumberedAntibody` objects; both CLI commands accept sequence
strings. No structure file or network connection is needed.

## Python

```python
import arpeggia

sequence = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGR"
    "FTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRGGYFDYWGQGTLVTVSS"
)
antibody = arpeggia.number_antibody(sequence, name="WT")
print(antibody)
print(antibody.cdr1, antibody.cdr2, antibody.cdr3)
for residue in antibody.residues:
    print(str(residue.position), residue.amino_acid, residue.input_index)

variant = arpeggia.number_antibody(sequence.replace("SSYAMS", "SSYALS"), name="Mutant")
alignment = arpeggia.align_antibodies([antibody, variant], reference_index=0)
print(alignment.format(width=100, color="never", rulers=False))
```

`NumberedAntibody` and its residue/reference records have read-only properties.
`input_sequence` retains the normalized complete input; `domain_span` is a
zero-based, half-open interval in that input. For example, `(5, 125)` identifies
input residues 6–125 in one-based notation: 120 supplied residues, excluding any
flanking tags or constant sequence. Imputation preserves this original interval.
`.sequence` contains the numbered residues only. `.fr1`, `.fr2`, `.fr3`, `.fr4`, `.cdr1`, `.cdr2`, and `.cdr3` are
strings under the selected CDR definition. Each residue carries its numbered
position, amino acid, original zero-based `input_index`, and region.

`chain` is `H`, `K`, or `L`. `confidence` is the numbering engine's normalized
profile-confidence heuristic, not a probability of correctness. `matched_profile_positions`
measures domain evidence without counting insertions. Diagnostics remain in
`.diagnostics`; Python also emits them as `UserWarning`.

## CLI

Set `sequence` to the amino-acid string above, then run:

```bash
arpeggia number-antibody "$sequence" --name WT --scheme imgt --species human,alpaca
arpeggia number-antibody "$sequence" --scheme martin --cdr-definition chothia --json
arpeggia align-antibodies "$sequence" "$sequence" --names WT,Replicate --reference-index 1
```

The positional arguments are sequence strings. `--names` assigns names in input
order; omitted or empty entries become `Seq001`, `Seq002`, and so on. Excess
names are an error. `--reference-index` is zero-based. `--json` emits plain
structured data, including tied references and imputation provenance.

## Numbering and CDR conventions

IMGT is the default. `scheme` / `--scheme` also accepts `chothia`, `martin`,
`aho`, and `kabat`. Chothia and Martin/enhanced Chothia use distinct rules,
including heavy FR3 insertion placement at H82 and H72, respectively.
`cdr_definition="auto"` uses the associated boundaries below. To override them,
explicitly provide both the numbering scheme and CDR definition.

| CDR definition | Heavy CDR1 / CDR2 / CDR3 | Light CDR1 / CDR2 / CDR3 |
| --- | --- | --- |
| IMGT | 27–38 / 56–65 / 105–117 | 27–38 / 56–65 / 105–117 |
| Martin/AbM | 26–35 / 50–58 / 95–102 | 24–34 / 50–56 / 89–97 |
| AHo structural loops | 25–40 / 58–77 / 109–137 | 25–40 / 58–77 / 109–137 |
| Kabat | 31–35 / 50–65 / 95–102 | 24–34 / 50–56 / 89–97 |
| Chothia (2021 consensus) | 26–32 / 52–56 / 96–101 | 26–32 / 50–52 / 91–96 |

Ranges are inclusive in each definition's native numbering, including insertions
at boundary positions. Mixed conventions transfer regions through the original
residue correspondence; they do not apply these numbers directly to labels from
another scheme. Thus `scheme="chothia"` uses Chothia numbering and consensus
regions, while `scheme="chothia", cdr_definition="martin"` retains Chothia labels
and assigns Martin/AbM regions.

Sources: [IMGT numbering](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html),
[Martin/AbM and Kabat definitions](https://pmc.ncbi.nlm.nih.gov/articles/PMC10939163/),
[AHo structural loops](https://pubs.rsc.org/en/content/articlehtml/2019/me/c9me00021f),
and the [upstream Chothia consensus table](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/numbering/chothia.rs#L16).

## Numbering without germline matching

For numbering and CDRs without V/J comparisons, use:

```python
numbered = arpeggia.number_antibody(sequence, match_germlines=False)
comparison = arpeggia.align_antibodies([numbered])
```

```bash
arpeggia number-antibody "$sequence" --no-germlines
arpeggia align-antibodies "$sequence" "$sequence" --no-germlines
```

Matching is enabled by default. Skipping it preserves numbering, regions,
confidence and domain detection; `species` has no effect in this mode. The
read-only `.germlines_searched` flag distinguishes skipped matching (`False`)
from a completed search (`True`), even if no reference qualified. Skipping sets
both `.v_match` and `.j_match` to `None` without a warning; other domain
diagnostics still apply.

The display says “Germline matching: skipped” and retains the input sequence,
rulers and CDR annotations, omitting the empty germline and operations rows.
Property access, formatting and JSON output never initiate matching.
`.impute()` raises `ValueError` if matching was skipped; number the input again
with matching enabled before imputing it. CLI `--no-germlines` conflicts with
`--impute`. A completed search with insufficient reference evidence retains the
usual imputation diagnostics. See the
[numbering-only benchmark](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md#optional-germline-matching)
for measured savings.

## Germline similarities

Offline references cover human, mouse, rat and rabbit H/K/L and alpaca/llama
heavy chains, including VHH references. Search uses all bundled species by default;
restrict it with `species="llama"`, `species=["rat", "rabbit"]`, or CLI
`--species rat,rabbit`. Both APIs accept these equivalent names, including the
same strains/subspecies. Quote Latin names in the shell: `--species "Homo sapiens"`.

| Common name | Latin alias |
| --- | --- |
| `human` | `Homo sapiens` |
| `mouse` | `Mus musculus` |
| `alpaca` | `Vicugna pacos` |
| `llama` | `Lama glama` |
| `rat` | `Rattus norvegicus` |
| `rabbit` | `Oryctolagus cuniculus` |

The llama supplement contains six V and five J references from IMGT protein
displays; it is a limited historical set, separate from the bulk snapshot.
The same generic H/K/L profiles number every input; the species option restricts
germline comparisons only. Adding species can change the best matches and ties
in an unrestricted search. It does not establish species-specific numbering accuracy.
Reported species identify references, not the organism of the input antibody.
See the [reference snapshot and attribution](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/README.md)
for source coverage, filtering, and known partial references.

`.v_match` and `.j_match` separately retain the highest qualifying local
BLOSUM62 score and all exact ties. V matching uses sequence through internal
IMGT 104; J matching uses sequence after it and requires FR4 evidence. A V hit
needs at least 50 known paired residues; J needs five known pairs within both
reference and input FR4. Insufficient evidence returns `None` with a diagnostic.

Each match's `.hits` groups references with identical sequence and coverage.
A hit stores its actual `SeqAlignment`, reference metadata, and these measures:

| Field | Meaning |
| --- | --- |
| `known_pairs`, `known_matches` | Paired canonical amino acids and identical pairs |
| `known_identity` | Known matches / known pairs |
| `reference_coverage` | Known pairs / known residues in the reference segment |
| `query_coverage` | Known pairs / known residues in the input segment |
| `query_input_start` | Input offset of the segment used by the alignment |
| `imgt_positions` | Reference residue positions; unavailable labels are `None` |

Ambiguous symbols and `U/O` do not count as known evidence. Metadata preserve
species, gene, allele, accession, and a stable snapshot-specific reference ID.
These are V/J similarities; they do not identify a unique ancestor or infer D.

## Terminal imputation

```python
partial = arpeggia.number_antibody(sequence[5:-3], name="Partial")
filled = partial.impute()
added = [residue for residue in filled.residues if residue.input_index is None]
```

Imputation returns a new object and fills only missing beginnings of FR1 or ends
of FR4. It preserves the supplied input, domain span, existing positions, internal
gaps, CDRs, and unknown input residues. Added residues have `input_index=None`
and `.imputed_from` reference IDs. Scores and coverage still describe the
supplied input. The read-only `.imputation_attempted` flag is `False` on fresh
numbering and `True` on the returned result, even when no residues were added.

Tied references must agree on a position's presence and known amino acid.
Conflicts or unavailable coverage remain unresolved with diagnostics. To choose
a specific tied reference, pass its exact `.references[i].id` to
`.impute(v_reference=..., j_reference=...)`. CLI equivalents are `--impute`,
`--v-reference`, and `--j-reference`; inspect `--json` to obtain IDs.

## Display and antibody alignments

The combined V/J germline is the reference above the supplied antibody query,
with the V gene on the left and J gene on the right. Each wrapped block has this order:

1. CDR marker
2. Germline ruler and reference sequence
3. Input ruler and query sequence
4. Operations relative to the germline

Only the representative V/J sequences are shown. The separate V/J summaries
identify them. `Additional tied V/J references` lists other gene/allele names
sharing the best score, grouped by species. It omits the displayed gene/allele
and collapses duplicate names; the record count covers only listed assignments.
The line is hidden when no additional names remain. All source records remain
in the result.

The numbered input defines the CDR bands across both rows. The top legend and
vertical bands use gray for CDR1, pink for CDR2 and cyan for CDR3. Bands cover markers, rulers and sequences, including gaps, but exclude
operations and name gutters. Yellow backgrounds mark imputed residues and the
`imputed residues` count. The count appears only after CLI `--impute` or Python
`.impute()`, including when the operation adds zero residues. Imputation
highlighting takes precedence over a CDR background at the same residue.

Rulers show one-based residue positions, with every tenth position right-aligned
above its residue. Endpoint numbers flank each row.
Original input coordinates survive imputation: prepending five inferred residues
leaves the first supplied residue at position 1. Imputed positions are blank.
The germline ruler counts continuously through V then J, ignoring alignment
gaps and the unknown junction. If V ends at 96, J starts at 97. This is a display
count; each stored V/J alignment retains its own source coordinates. The endpoint
follows the visible germline sequence with one space, then two spaces before the J name,
without padding to the input's right edge.
Canonical antibody numbering remains available through `.residues[i].position`.

Matches have blank operations; input insertions relative to the germline are
green `+`, input deletions red `-`, positive-BLOSUM62 substitutions blue `:`, and other
mismatches yellow `x`. Unknown V/J junctions have gray hyphens and blank operations.
Outer germline padding is blank where the input extends beyond germline coverage;
unmatched sequence tails are gray. These display choices preserve the local
V/J matches, scores, coverage and imputation evidence.

Python `.format(width=None, color="auto", rulers=True)` and CLI `--width`,
`--color auto|always|never`, and `--no-rulers` control the layout. Width includes
names and position gutters; automatic color follows terminal support and
`NO_COLOR`. Stored strings and JSON contain no ANSI escapes.

`AntibodyAlignment` stores the ordered union of `.positions`, corresponding
`.aligned_sequences`, original-order `.antibodies`, and `.reference_index`.
All antibodies must use one numbering scheme and be all heavy or all light;
K/L mixtures and differing per-row CDR definitions are allowed. Alignment follows
numbered positions and has no multiple-alignment score.

Python `NumberedPosition` objects compare and hash by `(number, insertion)`, so
residue positions can be looked up directly in `.positions`, sets and dictionaries.
The label carries no scheme or chain: retain that context across independent
results. Use alignment order rather than sorting position labels numerically.

The selected reference displays first, followed by the remaining inputs in their
original relative order. Germline rows are hidden, including for a one-antibody
alignment. `.format(reference_index=...)` changes one display without changing
stored rows, columns, or the default reference. The selected reference defines
the CDR bands across every sequence and ruler; the summary names that reference,
its numbering scheme and CDR definition. Each antibody retains its own region
annotations. Rulers follow each original input. If any included antibody has
undergone `.impute()`, the header `Total imputed residues: N` sums imputation
across every antibody, including the reference; otherwise the count is omitted.

## Scope and qualification

The adapter uses Immunum 1.3.1's Rust core. Inputs must be unaligned protein
strings of 30–10,000 residues, using the
[sequence alphabet](https://github.com/y1zhou/arpeggia/blob/master/docs/sequence-alignment.md#align-two-sequences).
Each must contain one recognizable
variable domain spanning FR1 through FR4; modest terminal FR1/FR4 truncations
are supported when the profile alignment spans IMGT 23–118, retaining coverage
of the FR1 and FR4 anchor positions. This preserves the framework context needed
for alignment and for both numbering and CDR conversion. Inputs cut past these
anchors are rejected; numbering cannot reliably infer the missing loop context.
Tags and constant-region tails are allowed but remain unnumbered.
Detected additional domains, severe partial domains, and insertions exceeding
the backend's single-letter representation produce errors. Recognition is a
heuristic and does not establish biological origin or numbering correctness.

Constant-region numbering, structure inputs, automatic multidomain handling,
multi-letter insertions, and species beyond the five bundled here are deferred. The
[qualification report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md)
records fixture agreement and limits; [ADR 0010](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md)
records the API decisions.
