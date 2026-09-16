# Separate sequence correspondence from RMSD evaluation

Extend [ADR 0008](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0008-cluster-structures-with-kabsch-and-k-medoids.md)
with optional sequence correspondence and rejection while preserving independent
fitting and evaluation populations. The [user guide](https://github.com/y1zhou/arpeggia/blob/master/docs/sequence-alignment.md)
defines arguments, result fields and display controls.

## Sequence alignment

`align_seqs()` accepts two unaligned amino-acid strings. Global alignment is the
default; local and query-full semi-global modes are selectable. Semi-global
consumes the entire second sequence with free first-sequence terminal overhangs.
FASTA parsing and general multiple sequence alignment are outside this API.

Use exact BLOSUM62 affine-gap optimization with matrix entries and gap costs
scaled together for integer scoring; report unscaled scores. Return one
deterministic optimum for the pinned backend. Structural RMSD never breaks
sequence-alignment ties.

Preserve original symbols for identity and edit distance when `U/O` score as
`C/K`. Report identity and coverage with both alignment-column and shorter-input
denominators so gaps and clipping remain interpretable. Edit distance compares
complete inputs independently of the chosen protein alignment. Similar
substitutions affect display, not identity or chemical equivalence.

The [alignment contract](https://github.com/y1zhou/arpeggia/blob/master/docs/sequence-alignment.md#align-two-sequences)
defines gap costs, alphabet validation, counts and empty-local-alignment results.

## Structural correspondence

Sequence correspondence applies to two-structure `rmsd`. Collection APIs still
require exact correspondence: pair-specific mappings must first be reconciled
with their shared coordinate layouts. Antibody numbering has a separate
[position-based contract](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0010-use-explicit-antibody-numbering-conventions.md).

Use observed chain sequences linked to coordinate residues. Declared residues
without coordinates cannot contribute atoms. Align complete chains before
applying residue ranges and atom selection. With alignment enabled, both
selectors address reference author numbering; otherwise selectors and exact
correspondence apply to both structures as in ADR 0008.

Explicit maps cover exactly the reference chains used by either selection,
assign unique query partners and disable inference. Otherwise score those
reference chains against all eligible query chains using shorter-against-longer
semi-global alignment. Maximize the summed raw score of a one-to-one assignment,
requiring every inferred pair score to be positive. A missing complete assignment
or tied optimum requires an explicit map. Unused query chains are allowed.
Identical reference scoring sequences, including `U/C` and `O/K` aliases, are
necessarily ambiguous and fail before pairwise scoring.

Final residue alignment uses the selected mode, with reference first and query
second, independently of chain-inference scoring. Report identity and coverage;
a positive score is not a universal homology threshold. Explicit maps can
override score-based inference, but refinement cannot guarantee removal of poor
sequence matches.

Nongap residue pairs include substitutions. Pair available backbone atoms across
substitutions; side chains require the same normalized amino-acid identity,
atom name and element. Scoring aliases do not make `U/O` chemically equivalent
to `C/K`. Report omitted atoms; exclude caps without sequence correspondence.

## Refinement and evaluation

`refine_cycles=0` performs one initial fit without rejection. Positive values
allow that many subsequent inspections, permanently rejecting atom pairs farther
than `refine_cutoff × current fitting RMSD` (default multiplier 2) and refitting
survivors. Stop when unchanged or effectively zero. If survivors cannot define
three non-collinear fitting pairs, fail with the surviving count.

Establish correspondence once. Refinement also works with exact correspondence
and does not infer structural matches as PyMOL `super` does. The final `rmsd`
evaluates every mapped pair selected by `rmsd_residues`, including rejected
fitting pairs, under the retained-pair transform. Evaluation never refits. This
keeps flexible-region deviations visible instead of reporting only a favorable
surviving-core score; `core_rmsd` reports that core separately.

## Results and presentation

Python `rmsd()` returns a read-only `RmsdResult` instead of a scalar. Record
initial/core/evaluation RMSDs, pair counts, inspection count, parameters and
chain/residue correspondence. Omit per-atom residuals and the transformation
to keep payloads compact. Rust, CLI JSON and Python share result structs, with
optional PyO3 annotations and binding methods rather than duplicate wrappers.

`SeqAlignment` retains normalized inputs, scoring settings, zero-based half-open
spans, statistics and equal-length gapped strings plus ASCII operations. These
recover residue-index pairs without a stored column vector. Names are metadata;
structural alignments append chain IDs. Stored data and JSON remain unstyled.

Share CLI/Python rendering while keeping stored data and JSON unstyled.
The [display contract](https://github.com/y1zhou/arpeggia/blob/master/docs/sequence-alignment.md#display-an-alignment)
owns names, operations, rulers and wrapping. RMSD's CLI includes chain alignments;
Python exposes them separately from its compact RMSD summary.

## Backend qualification

Hyalite 0.4.0 supplies exact modes and deterministic traceback without normal
dependencies. Arpeggia owns scoring semantics and result types. Its recent
introduction required independent qualification: 2,745 cases matched Biopython
scores, with traceback/span checks. Equal-score paths may differ across engines.

The [validation report](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/sequence-alignment.md)
records structural regressions, runtime and package-size measurements. Exact
RMSD retains its numerical behavior despite the changed result type.
[Sequence research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/protein-sequence-alignment-methods.md)
compares objectives and backends; the [PyMOL audit](https://github.com/y1zhou/arpeggia/blob/master/docs/research/pymol-superposition-and-rmsd.md)
establishes the relative rejection rule and core/full-evaluation distinction.

Backend sources: [Hyalite 0.4.0](https://docs.rs/crate/hyalite/0.4.0/source/)
and [Rust-Bio manifest](https://docs.rs/crate/bio/4.0.1/source/Cargo.toml.orig).
