# Use explicit antibody numbering conventions

Arpeggia will accept `chothia` as an alias for Martin/enhanced Chothia numbering,
as requested for the antibody-numbering API. This deliberately favors the
structurally corrected convention over compatibility with historical Chothia
outputs. The alias is an Arpeggia API choice: the original conventions remain
distinct in the [scheme authors' numbering service](https://www.bioinf.org.uk/abs/abnum/).

The initial CLI and Python APIs accept named amino-acid strings through
`number_antibody()` / `number-antibody`. A `NumberedAntibody` represents one
variable domain, retaining its input sequence and domain span. Terminal tags
and constant-region tails are allowed; additional detected variable domains
produce an error. Structure-file integration and constant-region numbering
are deferred. Initial partial-domain support is limited to truncation within
terminal FR1/FR4, retaining the intervening variable-domain core. More severe
truncations return an explicit error because length-based conversion can
renumber cut loops incorrectly. Supported partial domains retain coverage and
diagnostics; numbering itself does not fill missing sequence.
Ordered, read-only `.residues` records contain
the numbered position, amino acid, original input index and region. Derived
sequence properties expose `.sequence`, `.cdr1`/`.cdr2`/`.cdr3` and
`.fr1`/`.fr2`/`.fr3`/`.fr4`. Imputed residues have no original input index and
carry germline provenance. Chain, conventions and diagnostics remain available
on the object.

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

Display the top aligned V and J references together on one row, joined by
implicit gaps around CDR3. Show the V gene name on the left and J gene name on
the right to identify their distinct sources. The combined row is a display
of two reference matches, not an inferred ancestral sequence; their separate
scores, coverage, ties and provenance remain available.

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

An explicit imputation method returns a new `NumberedAntibody` and preserves
the original object and supplied sequence. It fills only missing chunks at the
beginning of FR1 or the end of FR4, using the closest aligned germline references.
It does not fill internal gaps or replace unknown input residues. A position
is filled only when tied references agree on its presence and amino acid;
otherwise it remains unresolved unless the caller selects a reference explicitly.
Record each inferred residue and its germline source so later alignments can
distinguish supplied and inferred residues.

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

The alignment CLI accepts positional sequence strings and comma-separated
`--names`, such as `--names WT,Mutant`. Omitted names become `Seq001`, `Seq002`,
and so on in input order. Missing or empty name entries use the default for
their sequence position; excess names are an error. Reuse JSON output, width,
color and ruler controls from sequence alignment.

Engine qualification compares RIOT, Immunum and AntPack on the
[AntPack test set](https://github.com/jlparkI/AntPack/tree/main/tests/test_data).
Pin the dataset and tool versions, retain failures, and distinguish agreement
with fixture labels from independent numbering accuracy. Report residue
numbering, domain coverage and failures separately.
The [comparison report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md)
records the completed 26,365-input run and unresolved qualification issues.

Use Immunum's native Rust core with default features disabled. Arpeggia owns
the result types, Python bindings and CLI/display formatting. Its four required
schemes and fixture agreement support this choice without establishing
independent accuracy. Initial recognition requires confidence at least 0.5
and 30 distinct matched profile positions, excluding query insertions.
This gate rejects trivial matches; it does not guarantee correct numbering.
Integration must correct the conversion/span defect and qualify long insertions and
multiple-domain detection. Pin the qualified core version and record any
upstream corrections. Reject unsupported insertion lengths before conversion
with a clear error; the initial integration does not extend Immunum's
single-letter insertion representation.

The V/J join's operation styling and the selected reference's display order
remain open.
Supporting evidence is in the
[numbering research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/antibody-numbering-schemes-and-tools.md).
