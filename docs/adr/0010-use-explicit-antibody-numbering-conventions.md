# Use explicit antibody numbering conventions

Use Immunum 1.3.1's native Rust core for variable-domain numbering; Arpeggia
owns results, germline matching, imputation, Python bindings and CLI rendering.
The [research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/antibody-numbering-schemes-and-tools.md)
compares alternatives and records upstream limitations. The
[qualification report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/antibody-numbering.md)
records fixture agreement and limits; the [user guide](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md)
defines arguments, fields and display details.

## Inputs and result objects

Accept named, unaligned amino-acid strings of 30–10,000 residues through
`number_antibody()` / `number-antibody`. Each `NumberedAntibody` represents
one variable domain and retains the normalized input, name and zero-based,
half-open domain span. Tags and constant tails remain unnumbered; additional
detected domains fail. Numbering does not fill missing sequence.

Recognition requires confidence ≥ 0.5 and 30 distinct matched profile positions,
excluding query insertions. Confidence is a heuristic, not a calibrated
probability. Require profile coverage spanning IMGT 23–118, including the
[FR1 Cys23 and FR4 W/F118 anchors](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html).
Loss of FR1 context can shift the raw alignment; this interval also contains
every count-based core conversion window in Immunum 1.3.1, including Kabat
heavy 24–40. AHo retains its upstream terminal-FR1 rule. Standalone V/J reference
segments are exempt from this whole-domain coverage gate.

Retain ordered, read-only residue records with typed position, amino acid,
original input index and region. Derive the numbered sequence and FR/CDR strings
from those records. Imputed residues have no original index and carry reference
provenance. Keep chain, conventions, coverage, confidence and diagnostics on the
result; retain private correspondence needed for rendering and imputation without
a second per-residue table. Position equality/hash compare the label, not its
scheme or chain context.

## Numbering and CDR definitions

IMGT is the default; Martin, AHo and Kabat are supported. Arpeggia deliberately
aliases `chothia` numbering to Martin/enhanced Chothia, favoring structural
corrections over historical Chothia output compatibility. The
[original conventions](https://www.bioinf.org.uk/abs/abnum/) remain distinct.

`cdr_definition="auto"` follows the numbering scheme. An explicit CDR override
requires an explicit scheme so mixed conventions are intentional. IMGT and Kabat
use their matching definitions; Martin uses AbM boundaries from the
[Martin group's study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10939163/).
AHo uses the published [structural-loop convention](https://pubs.rsc.org/en/content/articlehtml/2019/me/c9me00021f)
instead of Immunum's unverified table. Explicit Chothia uses Immunum's distinct
[2021 consensus definition](https://github.com/ENPICOM/immunum/blob/45bb70d34802cc592ebd86e685cc9f551885a2d6/src/numbering/chothia.rs#L16).
The [boundary table](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md#numbering-and-cdr-conventions)
and public docstrings cite these sources.

For mixed conventions, convert the same winning raw alignment into the CDR
definition's native scheme and transfer regions through input correspondence.
Do not apply another scheme's numeric boundaries directly or realign merely
to assign regions. Guard insertion capacities before conversion; do not silently
shorten an inconsistent domain or extend the backend's single-letter labels.
The adapter handles the qualified light-chain conversion/span defect explicitly.

## Germline matching

Bundle an attributed, versioned IMGT subset for human, mouse, rat and rabbit
H/K/L and alpaca/llama heavy chains/VHH. Llama uses a separately attributed
protein-display supplement because the bulk export omits it; retain its limited
coverage and unavailable J termini rather than treating it as a complete repertoire.
The same H/K/L profiles number all inputs;
species restrict reference matching only. Search all bundled species by default.
Report reference species, never presumed input origin or a unique ancestor.
Do not infer D segments. The [data record](https://github.com/y1zhou/arpeggia/blob/master/data/germlines/README.md)
defines filtering, coverage and licensing.

Compute separate V and J similarities once during numbering using local
BLOSUM62 with gap costs 10/0.5. V uses sequence through internal IMGT 104;
J uses sequence after it and requires FR4 evidence. Require 50 known paired V
residues or five known paired FR4 residues for J. Ambiguous symbols and `U/O`
do not count as known evidence. Retain every highest qualifying score tie and
its species/gene/allele/accession metadata, with identity and coverage separate
from score. Identical sequence-and-coverage references may share alignment work.
An unsupported segment returns no match and a diagnostic without losing numbering.

Rust/Python `match_germlines=false` / `False` and CLI `--no-germlines` skip matching
without changing numbering, CDRs or recognition; `species` is then unused.
`germlines_searched` distinguishes skipped searches from completed searches
without hits. Skipping leaves both matches absent, with no skip warning.
Reading, rendering and serialization never initiate matching. Imputation of an
unsearched result is an argument error; CLI `--no-germlines` conflicts with
`--impute`. Already-numbered objects with either search state can be aligned.

## Terminal imputation

Explicit `.impute()` returns a new antibody, preserving the input object,
sequence, domain span and supplied residue indices. Fill only missing beginnings
of FR1 and ends of FR4; internal gaps, CDRs and unknown input residues remain
unchanged. Reuse stored matches. Fill a position only when tied references agree
on its presence and known amino acid, unless the caller selects an exact tied
reference ID. Missing matches, unavailable endpoints and conflicts remain
unresolved with diagnostics. Every inferred residue retains its source IDs.
Retain `imputation_attempted`: zero added residues alone cannot distinguish an
unrequested operation from an attempted imputation that added nothing.

## Antibody alignments and display

`align_antibodies()` creates an ordered union of numbered positions and one
gapped sequence per antibody, retaining input row order and each row's CDR
annotations. Require one numbering scheme and all-heavy or all-light inputs;
kappa/lambda mixtures are allowed. Position correspondence is not a new MSA
calculation and has no fabricated alignment score.

A zero-based `reference_index` (default 0) controls comparison direction and
reference-defined CDR bands. A per-format override changes presentation, never
stored row order or columns. Show the reference first and hide germlines in
this view, including for one antibody. The CLI accepts positional strings and
comma-separated names, defaulting omitted/empty names by input position; excess
names and invalid reference indices fail.

For one numbered antibody, show the stitched V/J reference above the input query,
labeling both sources and listing additional tied gene/allele names in summaries.
Use germline-to-input operations directly: input insertions are + and deletions are -.
Keep every source record in the result. Unknown junctions
and missing outer coverage must remain distinguishable from alignment gaps.
Display padding changes never alter the local V/J matches or imputation evidence.

Input rulers preserve original coordinates and leave imputed residues blank.
The stitched germline ruler counts V then J continuously without gaps or junction
cells; stored V/J alignments retain separate source coordinates. This counter
does not imply a synthetic ancestor. Single-antibody CDR bands follow the numbered
input; multi-antibody bands follow the selected reference. Bands align across
rows, while imputation highlights retain each row's provenance. Shared sequence
rendering handles names, wrapping, color and rulers; JSON remains unstyled.
Exact layout and colors belong in the [display guide](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md#display-and-antibody-alignments).

## Qualification and deferred scope

Pin Immunum with default features disabled. Its Rust integration and strongest
fixture agreement across the four requested schemes support adoption; neither
fixture consensus nor successful numbering establishes independent accuracy.
The 26,365-input AntPack comparison retains rejections and distinguishes labels,
domain coverage and failures. Adapter qualification covers conversion boundaries,
recognition, partial domains, imputation and display, with runtime/package data
in the benchmark report.

Constant-region numbering, structure-file input, severe partial domains,
automatic multidomain handling, multi-letter insertions, additional species and
independent species-specific structural validation remain deferred.
