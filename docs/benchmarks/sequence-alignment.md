# Sequence alignment validation

## Backend qualification

Hyalite 0.4.0 was checked against Biopython 1.86 `PairwiseAligner` on
2026-09-09. Both used BLOSUM62 and affine costs `open + (length - 1) * extend`.
Hyalite matrix entries and costs were multiplied by 100; reported scores were
unscaled. Semi-global alignment consumed the whole second sequence while allowing
free first-sequence terminal overhangs.

The 2,745 cases combined five curated pairs and 300 pseudorandom pairs
(seed 902, lengths 1–65, alphabet `ARNDCQEGHILKMFPSTWYVBZX`), three modes,
and opening/extension costs `(10, 0.5)`, `(1.25, 0.25)`, and `(0.5, 0.5)`.
Every optimal score matched the reference. Each traceback reconstructed its
score and consumed its declared spans; global and semi-global endpoints obeyed
the requested policy. Equal-score paths need not match across implementations.

The release-mode Hyalite subprocess processed these cases in 0.426 s on the
local Linux x86-64 host. This includes input/output and scoring-object creation;
it is a qualification-run observation, not a comparative speed benchmark.
Full traceback is scalar; SIMD score-only claims do not describe this workload.

Pinned Hyalite adds one package to the lockfile and has no normal dependencies.
Serde and serde_json, now direct dependencies for result serialization, were
already present transitively. Final wheel measurements will include the complete
feature rather than attributing all growth to Hyalite.

Repository regressions cover directional semi-global spans, gap costs and counts,
empty local alignments, original-symbol identity for U/O scoring aliases, invalid
inputs, and identical full-matrix/checkpoint traceback on repeated sequences.
The qualification harness and its Python environment are not production dependencies.
