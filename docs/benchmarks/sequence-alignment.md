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
already present transitively.

Repository regressions cover directional semi-global spans, gap costs and counts,
empty local alignments, original-symbol identity for U/O scoring aliases, invalid
inputs, and identical full-matrix/checkpoint traceback on repeated sequences.
The qualification harness and its Python environment are not production dependencies.

## Structural correspondence regressions

A 4,096-matrix exhaustive check validates the assignment solver against all
eligible two-chain assignments with negative, zero, and positive scores.
Fixtures cover renumbering, substitutions, directional reference selections,
unused query chains, ambiguous homomers, explicit overrides, gaps, and chemical
side-chain exclusions. Existing exact-correspondence numerical tests still pass.

The symmetric 20-pair refinement fixture from the PyMOL research retains 16
fitting pairs at core RMSD 0.1 after two rejection passes and one unchanged
inspection. Evaluation still includes all 20 pairs at RMSD approximately 1.925617.
A threshold that removes every fitting pair fails with the surviving count.
An unchanged refinement inspection preserves the 0.671302 Å RMSD of a reflected
four-point fixture; impossible reference/query chain counts fail before scoring.

## Python performance and package size

Release wheels for baseline `e3a1677` (v0.9.2) and measured feature build `aded22f` used the same
Rust 1.96 toolchain and CPython 3.13 on Linux x86-64 (Ryzen 9 9950X3D).
The table reports the median of five fresh-process medians, with three warm-up
calls per process and alternating case order.

| Operation | Median time (ms) |
| --- | ---: |
| Baseline exact RMSD | 0.959 |
| Feature exact RMSD | 0.924 |
| Feature sequence-aligned RMSD | 0.945 |
| Global sequence alignment, 80 residues | 0.023 |
| Global sequence alignment, 256 residues | 0.279 |
| Global sequence alignment, 1,024 residues | 4.609 |

RMSD compares the repository's 1UBQ fixture with itself, using default selections
and 100 timed calls per process. Timings include parsing and Python result
construction. All results were numerically zero. The small differences show no
meaningful regression on this input; they do not establish a speedup or predict
large multichain performance.

Sequence inputs repeat `ACDEFGHIKLMNPQRSTVWY`, truncated to the stated length;
the second input inserts `GG` at the midpoint. Global defaults apply, with
125, 39, and 9 timed calls per process respectively. Each result has edit
distance two. Timings include alignment, metrics, and Python result handling.

| Artifact | Baseline bytes | Feature bytes | Growth |
| --- | ---: | ---: | ---: |
| Compressed wheel | 9,487,996 | 9,596,012 | 1.14% |
| Native extension | 34,653,584 | 34,910,936 | 0.74% |

Growth includes the complete feature, bindings, and serialization. Validation
passed 205 Rust tests (including CLI and doctests), 17 Python tests against the
fresh wheel with warnings treated as errors, Python type checking, and all
pre-commit checks.

## Gapped results and display validation

The artifact measurements above cover the alignment engine and API before display
rendering and expanded help; they are not final release-size measurements.
The display-stage build passed 211 Rust tests and 18 Python tests against the rebuilt editable
extension, with Python warnings treated as errors. Regressions cover gapped
strings and recovered residue indices, insertion/deletion orientation,
positive-score substitution markers without changing identity, independent
rulers across gaps, named Unicode labels, narrow wrapping, gray clipping, empty local results, and plain
JSON. CLI and Python were also checked in a 40-column pseudo-terminal: both
detected width, enabled color automatically, and respected `NO_COLOR`.
