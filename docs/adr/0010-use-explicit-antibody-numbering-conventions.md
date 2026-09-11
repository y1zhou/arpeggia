# Use explicit antibody numbering conventions

Implemented and qualified on 11 September 2026. The
[implementation record](https://github.com/y1zhou/arpeggia/blob/master/docs/research/antibody-numbering-schemes-and-tools.md#10-implementation-plan)
records milestones; the [qualification report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md#arpeggia-adapter-qualification)
records validation and limits.

## Inputs and result objects

Arpeggia accepts `chothia` as an alias for Martin/enhanced Chothia numbering,
as requested for the antibody-numbering API. This deliberately favors the
structurally corrected convention over compatibility with historical Chothia
outputs. The alias is an Arpeggia API choice: the original conventions remain
distinct in the [scheme authors' numbering service](https://www.bioinf.org.uk/abs/abnum/).
IMGT is the default; Martin, AHo and Kabat are also supported.

The initial CLI and Python APIs accept named amino-acid strings through
`number_antibody()` / `number-antibody`. A `NumberedAntibody` represents one
variable domain, retaining its input sequence, name and zero-based, half-open
domain span. Terminal tags and constant-region tails are allowed; additional
detected variable domains produce an error. Structure-file integration and constant-region numbering
are deferred. Initial partial-domain support is limited to truncation within
terminal FR1/FR4, retaining the intervening variable-domain core. More severe
truncations return an explicit error because length-based conversion can
renumber cut loops incorrectly. Supported partial domains retain coverage and
diagnostics; numbering itself does not fill missing sequence.
Ordered, read-only `.residues` records contain
the numbered position, amino acid, original input index and region. Derived
sequence properties expose `.sequence`, `.cdr1`/`.cdr2`/`.cdr3` and
`.fr1`/`.fr2`/`.fr3`/`.fr4`. Imputed residues have no original input index and
carry germline provenance. Chain, conventions, confidence, matched-profile
coverage and diagnostics remain available on the object. Private correspondence
and reference-coverage data may support rendering and imputation; derived views
do not require duplicate per-residue tables.

## Germline matching

The numbering call computes germline matches once, returning separate V and J
reference similarities, coverage and tied references for display and imputation.
It does not reconstruct a unique ancestral antibody or
infer a D segment. References ship offline as a versioned, attributed IMGT
subset covering human and mouse H/K/L and alpaca heavy-chain/VHH references.
Search covers all bundled species by default and accepts an explicit species
restriction. Report matched-reference species rather than presumed input origin.

Rank V and J references separately by local BLOSUM62 alignment score, using
the sequence module's gap costs of 10/0.5. The V search uses observed sequence
through internal IMGT position 104; J searches after it and must reach FR4.
Retain exact score ties and report identity and coverage separately. Initial
qualification gates require 50 known paired V residues and five known paired
FR4 residues for J. Ambiguous residues do not count as known evidence. If a
segment lacks enough evidence, retain the numbered antibody with `v_match=None`
or `j_match=None` and a diagnostic; imputation uses only supported references.
Preserve species, gene, allele and accession metadata for all tied references.
Choose a deterministic display representative without discarding ties. Identical
reference sequences with identical coverage may share alignment work.

## Numbered-antibody display

Display the top aligned V and J references together on one row, joined by
gray hyphens around the uncovered V/J junction, with blank operation markers.
These cells represent unavailable reference sequence and must remain distinct
from aligned insertions or deletions. Show the V gene name on the left and
J gene name on the right to identify their distinct sources. The combined row
displays two reference matches, not an inferred ancestral sequence; their separate
scores, coverage, ties and provenance remain available.

Use three distinct CDR background colors with regions identifiable in plain
output too. Reuse sequence-alignment foreground conventions: blank for a match,
green `+` for insertion, red `-` for deletion, blue `:` for a positive-BLOSUM62
substitution, and yellow `x` for other mismatches. Unmatched tails are gray.
Stored sequences and JSON contain no ANSI styling.

CLI and Python displays infer terminal width and color support; explicit width,
color and ruler controls follow `SeqAlignment`. Python supports
`.format(width=None, color="auto", rulers=True)`. The width includes both name
gutters and position labels. Rulers follow canonical numbered positions;
imputed residues and stitched references do not create fictitious input offsets.
Reuse the existing renderer through a small shared layout helper, without
constructing a `SeqAlignment` whose optimal-alignment score would be misleading.

## Numbering and CDR definitions

`cdr_definition="auto"` selects the CDR convention associated with the chosen
numbering scheme. A caller choosing a different CDR definition must explicitly
provide both arguments. This keeps the usual scheme/region pairing convenient
while making a mixed convention intentional. IMGT and Kabat use their matching
definitions. Martin uses AbM/Martin boundaries, as in the
[Martin group's study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10939163/).
AHo uses CDR1 25–40, CDR2 58–77 and CDR3 109–137 from the published
[structural-loop convention](https://pubs.rsc.org/en/content/articlehtml/2019/me/c9me00021f).
Public docstrings must cite these definitions. Explicit
`cdr_definition="chothia"` selects the distinct 2021 consensus Chothia definition
used by Immunum: heavy 26–32 / 52–56 / 96–101 and light 26–32 / 50–52 / 91–96.
Thus `scheme="chothia", cdr_definition="auto"` resolves to Martin/AbM, while an
explicit Chothia CDR override retains Martin numbering and uses Chothia regions.
See the [cited upstream table](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/numbering/chothia.rs#L16).
For mixed conventions, convert the same winning raw alignment to the CDR
definition's native scheme and transfer regions through input correspondence.
Do not apply one scheme's numeric boundaries to another scheme's labels or
repeat sequence alignment solely to assign regions.

## Terminal imputation

An explicit `.impute()` method returns a new `NumberedAntibody` and preserves
the original object and supplied sequence. It fills only missing chunks at the
beginning of FR1 or the end of FR4, using the closest aligned germline references.
It does not fill internal gaps or replace unknown input residues. A position
is filled only when tied references agree on its presence and known amino acid;
otherwise it remains unresolved unless the caller selects a reference explicitly.
An unavailable match or uncovered reference end leaves that end unresolved with
a diagnostic. Ambiguous reference residues cannot supply imputed residues.
Record each inferred residue and its germline source so later alignments can
distinguish supplied and inferred residues.

## Antibody alignments

`align_antibodies()` / `align-antibodies` produces an `AntibodyAlignment` from
one or more numbered antibodies in one numbering scheme. Rows must all be heavy
chains or all be light chains; kappa/lambda mixtures are allowed. Reject
incompatible inputs explicitly. Each row retains its CDR definition because
column correspondence follows numbered positions rather than region labels.
Expose `.antibodies`, the ordered union of `.positions`, and corresponding
`.aligned_sequences`, preserving input row order. Use a zero-based
`reference_index`, defaulting to 0, to identify the antibody used for display
comparisons. Accept it in `align_antibodies()` and as `--reference-index` in the
CLI; `.format(reference_index=...)` can override it for one rendering. Changing
the reference changes comparison direction and colors without renumbering or
altering stored columns and rows. Reject out-of-range indices.
Display the selected reference first, followed by the other antibodies in their
original relative order. Hide germline sequences in this view, including a
single-antibody alignment. The reference retains its CDR backgrounds but has
neutral foreground coloring and no self-comparison operation row. Other rows
show operations relative to it. These presentation choices do not reorder
`.antibodies` or `.aligned_sequences` or introduce a multiple-alignment score.

The alignment CLI accepts positional sequence strings and comma-separated
`--names`, such as `--names WT,Mutant`. Omitted names become `Seq001`, `Seq002`,
and so on in input order. Missing or empty name entries use the default for
their sequence position; excess names are an error. Reuse JSON output, width,
color and ruler controls from sequence alignment.

## Backend qualification

Engine qualification compares RIOT, Immunum and AntPack on the
[AntPack test set](https://github.com/jlparkI/AntPack/tree/main/tests/test_data).
Pin the dataset and tool versions, retain failures, and distinguish agreement
with fixture labels from independent numbering accuracy. Report residue
numbering, domain coverage and failures separately.
The [comparison report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md)
records the 26,365-input engine comparison and qualified adapter results.

Use Immunum 1.3.1's native Rust core with default features disabled. Arpeggia owns
the `ab_numbering` module, result types, Python bindings and CLI/display
formatting. Its four required schemes and fixture agreement support this choice
without establishing independent accuracy. Initial recognition requires confidence at least 0.5
and 30 distinct matched profile positions, excluding query insertions.
This gate rejects trivial matches; it does not guarantee correct numbering.
The adapter corrects the conversion/span defect and guards long insertions and
multiple-domain detection. Pin the qualified core version and record any
upstream corrections. Reject unsupported insertion lengths before conversion
with a clear error; the initial integration does not extend Immunum's
single-letter insertion representation.

The [numbering research and implementation plan](https://github.com/y1zhou/arpeggia/blob/master/docs/research/antibody-numbering-schemes-and-tools.md#10-implementation-plan)
record the integration sequence, validation criteria and remaining engineering
choices. Constant-region numbering, structure-file integration and additional
species beyond human, mouse and alpaca remain deferred.
