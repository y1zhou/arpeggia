# Separate sequence correspondence from RMSD evaluation

This decision extends [ADR 0008](0008-cluster-structures-with-kabsch-and-k-medoids.md)
with optional sequence correspondence while preserving independent fitting and
evaluation populations.

## Sequence alignment

`align_seqs()` returns a `SeqAlignment` for two unaligned amino-acid strings
for Rust, Python, and CLI use. Global alignment is the default; local and
semi-global modes are selectable. Semi-global consumes the entire second
sequence and permits free terminal overhangs of the first. Alignment uses exact
affine-gap optimization with BLOSUM62. Configurable positive gap costs default
to opening 10 and extension 0.5, with cost `open + (length - 1) * extend`.
Costs have at most two decimal places and opening cost must be at least extension
cost. Unsupported precision is rejected rather than rounded. Matrix entries and
gap costs are scaled consistently for exact integer scoring; reported scores are
unscaled. Nonfinite costs and unrepresentable score ranges fail explicitly.
FASTA parsing and multiple sequence alignment are outside this feature.

Identity and paired-residue coverage each have two normalizations. For `M`
identical pairs, `P` nongap pairs, `A` alignment columns including gaps, and `S`
the shorter full input length, identity is `M/A` and `M/S`; coverage is `P/A`
and `P/S`. Retain the underlying counts. Unaligned terminal segments do not count
toward `A`; substitutions count toward `P` but not `M`. Empty local alignments
have undefined alignment-length ratios and zero shorter-input ratios.

Gap-residue and gap-run counts are distinct. Edit distance is full-input
Levenshtein distance, independent of the protein-scored alignment; local clipping
does not shorten its inputs.

Lowercase is normalized. Standard amino acids and `B/Z/X` are accepted. `U/O`
remain distinct input symbols for identity but score as `C/K`, with a diagnostic.
Input gaps, stop symbols, unsupported characters, and empty strings are rejected.
A local alignment without a positive-scoring match returns an empty alignment:
score zero, the empty-alignment ratios defined above, and full-input edit distance. One deterministic optimal traceback is returned for the pinned
backend; structural RMSD is not used to break sequence-alignment ties.

## Structural correspondence

Sequence correspondence applies to two-structure `rmsd`. Extending pairwise-table
and clustering APIs requires further work: pair-specific correspondence must first
be reconciled with their shared coordinate layouts and comparison semantics.
Future antibody numbering can consume sequence correspondence but requires its
own domain and numbering definitions.

With `--align-seqs`, structural mapping uses observed sequences, retaining their
links to coordinate residues. Declared residues without coordinates cannot
participate in superposition. Single-chain inputs are paired automatically;
explicit chain mappings are supported. For multi-chain inputs without identical
homomer chains, infer mappings from all-to-all semi-global alignments, consuming
the shorter chain against the longer. Maximize summed raw scores over a
one-to-one chain assignment. A tied optimum requires explicit mapping. Every
reference chain used by either residue selection requires a partner; unused
query chains are allowed. Inferred pair scores must be positive, with identity
and coverage reported; no universal homology threshold is imposed. Explicit
mapping can override score-based inference. Infer only for reference chains used
by either selection, against all eligible query chains. Explicit maps cover all
relevant reference chains, use unique query partners, and disable inference.

Chain inference and final residue alignment are separate: inference uses
shorter-against-longer semi-global scores, while final alignment uses the selected
mode with reference first and query second. In alignment mode, both residue
selectors address reference author numbering and map to query residues after
complete observed chains have been aligned. Without alignment, existing exact
correspondence and selector behavior remain.

Nongap residue pairs include substitutions. Available backbone atoms pair across
substitutions; side-chain pairing requires the same normalized amino-acid
identity, followed by matching atom names and elements. The `U/O` scoring aliases
do not establish chemical equivalence with `C/K`. Omitted
atoms are reported, and caps lacking sequence correspondence are excluded.
Geometric refinement cannot establish chemical equivalence or guarantee rejection
of a poor sequence match.

## Refinement and evaluation

Refinement is independent of sequence alignment and can use exact atom
correspondence. It does not infer structural correspondence as PyMOL `super` does.
`refine_cycles=0` preserves one initial fit without rejection. A positive value
allows that many subsequent rejection/refit cycles. Rejection is permanent and
atom-wise, using a configurable multiplier of the current fitting RMSD, initially
2. Sequence correspondence is established once, before refinement. Stop early
when no pairs are rejected or fitting RMSD is effectively zero. If rejection
leaves fewer than three non-collinear fitting pairs, fail with the surviving count
rather than silently falling back to the preceding fit.

The final RMSD evaluates all mapped pairs selected by `rmsd_residues`, including
pairs rejected from fitting, under the final retained-pair transform. Evaluation
never triggers another fit. This preserves flexible-region deviations rather
than reporting only a favorable surviving-core RMSD.

## Result contract

`rmsd()` returns a `RmsdResult` instead of a scalar, accepting the breaking API
change. Its `rmsd` is final evaluation RMSD; `core_rmsd` is the RMSD of retained
fitting pairs after rejection. It also records initial fitting RMSD, initial and
retained fitting counts, evaluation count, performed cycles, and chain/residue
correspondence. These populations can differ when fit and evaluation selections
differ.

Both result types expose read-only Python attributes; Python continues emitting
scientific warnings. `SeqAlignment` retains normalized inputs, scoring settings,
aligned spans, residue mapping, and statistics. Per-atom residuals, atom-pair
records, and the final transformation are omitted to keep results compact.

The CLI returns a concise detailed summary by default, with explicit JSON output
for parameters, mappings, and statistics.

### Gapped sequences and display

Sequence names are display metadata, separate from sequence data and residue
correspondence. `reference_name` and `query_name` default to "Reference" and
"Query"; structural alignments append their chain IDs. Names are included in
results and JSON, with label padding based on visible terminal width.

`SeqAlignment` exposes equal-length `aligned_reference`,
`aligned_query`, and `operations` strings. Gapped strings make downstream use
direct; original inputs and zero-based, half-open spans retain enough information
to reconstruct index pairs internally. Strings contain only the scored alignment,
using `-` for gaps and no color escapes. Operations describe reference-to-query
changes: space for match, `+` for insertion, `-` for deletion, and `x` for
substitution.

Render reference and query rows followed by an unlabeled operation row. Insertions are green,
deletions red, and mismatches yellow, coloring both sequence cells and the marker;
matches are uncolored. Clipped tails are gray with blank operation cells, with
prefixes right-aligned against the scored alignment and suffixes left-aligned
after it. Global terminal gaps remain scored operations. Empty local results
show both inputs gray and state that no positive-scoring alignment exists.

Wrap complete output into blocks within terminal width, falling back to 80
columns, or an explicit width including labels and positions. Reject widths
that cannot fit labels and one residue. CLI and Python object displays enable
color automatically when the terminal supports it; stored fields and JSON remain
plain. Explicit color controls and `NO_COLOR` support remain available.

Display positions are one-based input coordinates, counting residues but not gaps
or padding. Each sequence row shows its start and end positions; rulers mark
every tenth residue with the last digit aligned to that residue's column.
Rows without residues omit endpoint numbers. Structural sequence displays use
observed-sequence positions; author numbering remains in the residue mapping.
Each ruler sits immediately above its sequence row, with operations below both.
CLI `--no-rulers` and Python `.format(rulers=False)` hide rulers while retaining
start/end positions; rulers are shown by default.
Python `repr()`, `str()`, and `.format(width=None, color="auto", rulers=True)` use automatic
color. Explicit `color="never"` produces plain text; `color="always"` overrides
terminal detection and `NO_COLOR`.

The RMSD CLI shows each chain alignment. Python `RmsdResult` remains compact,
with full displays available through its individual chain alignment objects.

## Backend qualification

Hyalite 0.4.0 provides exact alignment modes and deterministic traceback without
normal dependencies. Arpeggia owns the public result types and scoring semantics.
Its recent introduction warranted independent qualification before adoption:
2,745 cases matched Biopython scores, with tracebacks checked for score and span
consistency. Tied optima require deterministic results within the pinned backend,
not identical paths across implementations.

Structural regressions cover renumbering, chain assignment and ambiguity,
substitutions, missing atoms, independent fit/evaluation selections, refinement
failure, and evaluation of pairs rejected from fitting. Exact-correspondence RMSD
retains its numerical behavior despite the new result type. The
[validation report](../benchmarks/sequence-alignment.md) records the checks and
measured runtime and package sizes.

Backend evidence: [Hyalite 0.4.0 source](https://docs.rs/crate/hyalite/0.4.0/source/)
and [Rust-Bio manifest](https://docs.rs/crate/bio/4.0.1/source/Cargo.toml.orig).

The [sequence research](../research/protein-sequence-alignment-methods.md)
explains scoring and end-gap objectives; the
[PyMOL audit](../research/pymol-superposition-and-rmsd.md) documents the relative
rejection rule and the distinction between core and full-pair evaluation.
