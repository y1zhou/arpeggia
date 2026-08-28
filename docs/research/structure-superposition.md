# Protein structure superposition and RMSD

Research date: 2026-08-28

## Decision summary

Kabsch/SVD and quaternion characteristic polynomial (QCP) solve the same
least-squares rigid-superposition problem. Given the same paired atoms,
weights, centering convention, and proper-rotation constraint, they should
produce the same minimum RMSD up to floating-point error. The solver should
therefore remain an implementation detail rather than a user-visible
scientific-method option.

The accepted solver direction after the design interview is:

1. Define and test one small, serial `f64` Kabsch kernel whose inputs are two
   already-paired `(n, 3)` coordinate arrays with equal shapes.
2. Use the existing `nalgebra` dependency for the 3-by-3 SVD and enforce a
   proper rotation. The public result is RMSD; the transform remains internal.
3. Parallelize independent structure pairs when constructing an ensemble
   matrix, not atoms within one fit.
4. Do not implement QCP or add a third-party superposition crate in this
   version. The comparison below records why QCP remains a credible future
   optimization rather than a second public method.

## Scientific equivalence

For paired, centered coordinate matrices, both methods minimize the weighted
sum of squared distances under a rotation. Kabsch obtains the optimal rotation
from an SVD of the 3-by-3 cross-covariance matrix. QCP expresses the same
objective through a 4-by-4 quaternion key matrix, locates its largest eigenvalue
from the characteristic polynomial, and obtains RMSD from that eigenvalue.
The 2010 QCP extension obtains the corresponding rotation quaternion from the
adjoint matrix. The QCP derivation explicitly includes a diagonal atom-weight
matrix. See [Kabsch 1976](https://doi.org/10.1107/S0567739476001873),
[Kabsch 1978](https://doi.org/10.1107/S0567739478001680),
[Theobald 2005](https://doi.org/10.1107/S0108767305015266), and
[Liu, Agrafiotis, and Theobald 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2958452/).

Neither solver should permit reflection for protein superposition. A reflection
can reduce a least-squares residual by reversing chirality, but it is not a
physical rigid-body rotation. Kabsch must explicitly correct the SVD result so
that `det(rotation) = +1`; a unit quaternion represents a proper rotation by
construction. Umeyama's analysis is a useful warning that an unchecked SVD
solution can return a reflection for corrupted data
([Umeyama 1991](https://web.stanford.edu/class/cs273/refs/umeyama.pdf)).

Scaling must also remain disabled. Similarity-transform variants of Umeyama
estimate scale, but changing molecular bond lengths is not part of structural
superposition and would change the meaning of RMSD.

## Comparison

| Property | Kabsch/SVD | QCP |
| --- | --- | --- |
| Scientific objective | Minimum least-squares RMSD under a proper rigid rotation | Same |
| Main fixed-size solve | SVD of a 3-by-3 cross-covariance matrix | Largest root of a quartic; optionally recover a quaternion rotation |
| RMSD without transform | Still normally computes the SVD | Avoids decomposition and can return RMSD directly |
| Transform | Directly available from the SVD factors | Available through the 2010 adjoint-matrix extension |
| Atom weights | Incorporated in centroids and covariance | Explicitly supported by the published weighted formulation |
| Reflection handling | Requires determinant correction | Unit quaternion yields a proper rotation |
| Implementation fit here | `nalgebra` 0.35 is already a dependency | No suitable focused Rust dependency found; requires a carefully validated local implementation |
| Extension beyond 3D | SVD formulation generalizes naturally | Quaternion formulation is specifically three-dimensional |
| Typical ecosystem | General crystallography, geometry, and structural toolkits | High-throughput molecular trajectory and RMSD workloads |

## Performance evidence and its limit

The 2010 paper reports that QCP recovered optimal rotations about 20 times faster
than a 4-by-4 Householder/QL eigen decomposition and produced the same rotations
within floating-point error. However, its timing explicitly excludes construction
of the coordinate inner-product matrix because every method needs that work
([paper text and Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC2958452/)).
The headline factor is therefore not an end-to-end protein-pair speedup.

For `S` structures with `N` selected atoms, an all-pairs matrix still requires
`O(S^2 N)` coordinate products with either solver, plus `O(S^2)` fixed-size
solves. Centered coordinates, centroids, and each structure's self inner product
can be cached once, but the cross inner product cannot. Consequences:

- QCP should matter most for C-alpha or other small selections, where the
  fixed-size solve is a larger fraction of each comparison.
- For full-atom multimers, the `O(N)` cross-covariance pass may dominate, making
  the difference between a 3-by-3 SVD and QCP modest.
- Pair-level Rayon parallelism is likely more consequential than solver choice,
  but nested parallelism should be avoided.
- A benchmark must include atom access, cached centering, covariance accumulation,
  the solve, and RMSD production. A microbenchmark of only SVD versus QCP would
  repeat the limitation of the published timing.

QCP is nevertheless well established in performance-oriented molecular tools.
MDTraj uses Theobald QCP for its default superposed RMSD and parallelizes across
frames ([MDTraj RMSD API](https://mdtraj.readthedocs.io/en/latest/api/generated/mdtraj.rmsd.html));
MDAnalysis exposes weighted `float32`/`float64` QCP RMSD and rotation while doing
the QCP arithmetic internally in double precision
([MDAnalysis QCP API](https://docs.mdanalysis.org/stable/documentation_pages/lib/qcprot.html)).

## Parallelization boundary

For an ensemble of `S` structures there are `S(S-1)/2` independent pairs. At
only 1,000 structures this is already 499,500 tasks, so pair-level parallelism
provides ample work for a fixed Rayon pool. Each worker should run the serial
Kabsch kernel for its assigned pairs and write to disjoint locations in the
packed matrix. Contiguous matrix chunks can preserve useful row locality while
avoiding an auxiliary vector of pair indices.

Parallelizing atoms within each pair is a poorer initial boundary. Every pair
would repeatedly schedule a reduction for a small 3-by-3 covariance, followed
by a constant-size SVD that has no useful internal parallel work. Nesting that
inside a parallel pair loop adds scheduling and reduction overhead, makes
floating-point accumulation order depend on scheduling, and risks using an
ambient Rayon pool differently from the explicitly configured Arpeggia pool.
Rayon recommends sequential inner iteration when an outer parallel iterator
already supplies the concurrency
([`flat_map_iter` rationale](https://github.com/rayon-rs/rayon/blob/main/RELEASES.md)).

The implementation should therefore install one existing Arpeggia-configured
Rayon pool around packed-matrix construction, split the output into disjoint
chunks, and calculate each pair serially in canonical atom order. This keeps
each RMSD numerically deterministic regardless of thread scheduling. The
standalone two-structure `rmsd` operation remains serial. Per-atom threading is
a future optimization only if benchmarks on exceptionally large selections
show that the scalar operation needs it.

Structure preparation is a separate bounded parallel phase. Parse the first
structure serially to establish the reference atom identities and selected atom
count, perform the second memory check, then prepare the remaining structures
with `min(num_threads, 8)` workers. Collect results in canonical input order so
warnings and failures remain deterministic. Drop each full parsed structure
after retaining only its centered selected coordinates. The eight-worker cap is
an I/O-pressure heuristic and must be benchmarked rather than presented as a
hardware guarantee.

## Numerical robustness

Use `f64` internally for both methods. PDB and mmCIF coordinate precision does
not justify accepting extra cancellation in an all-pairs distance matrix, and
MDAnalysis likewise keeps the QCP calculation in double precision even when
input coordinates are `float32`
([MDAnalysis QCP API](https://docs.mdanalysis.org/stable/documentation_pages/lib/qcprot.html)).

Kabsch has the advantage of relying on `nalgebra`'s maintained SVD. Its fallible
`try_new` API permits an explicit convergence limit rather than an unbounded
solve ([nalgebra SVD API](https://docs.rs/nalgebra/0.35.0/nalgebra/linalg/struct.SVD.html)).
QCP avoids a general decomposition, and the authors report rapid, stable
Newton-Raphson convergence from the self-inner-product upper bound. The 2010
paper reports roughly five iterations for relative precision `1e-6` and more
than one billion tested fragment superpositions
([Liu et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2958452/)).

QCP still needs explicit failure policy and adversarial tests. Biopython's current
implementation retains a bounded Newton iteration and fallbacks when candidate
adjoint columns are too small
([Biopython QCP source](https://github.com/biopython/biopython/blob/master/Bio/PDB/qcprot.py));
MDAnalysis has previously fixed a case where its QCP routine returned no RMSD
([MDAnalysis changelog](https://github.com/MDAnalysis/mdanalysis/blob/develop/package/CHANGELOG)).
Those are reasons to validate a local implementation against Kabsch rather than
to treat the formula as automatically infallible.

Required numerical tests for the production Kabsch solver:

- identical coordinates and pure translations/rotations;
- noisy coordinates with a known transform;
- mirrored coordinates, verifying `det(rotation) = +1`;
- very large coordinate offsets, verifying that centering avoids loss of
  precision;
- repeated, collinear, and coplanar points, including the accepted degeneracy
  failures;
- near-zero RMSD, where an algebraic residual can become slightly negative
  from rounding;
- the accepted uniform-weight convention;
- symmetry: `rmsd(A, B)` approximately equals `rmsd(B, A)`;
- agreement across supported thread counts and randomized rigid transforms
  within a stated absolute tolerance.

A transform is non-unique for degenerate point sets even when the minimum RMSD
is well defined. The accepted public behavior rejects empty selections and
requires at least three non-collinear points, giving an unambiguous protein
superposition boundary.

## Atom pairing and subsets

Chain, residue-range, and atom-type selection happens before the numerical
solver. It must produce two coordinate arrays in one deterministic biological
identity order, not merely in each file's record order. For the proposed first
version, a safe identity key includes at least model, chain ID, residue serial,
insertion code, residue name, atom name, and the selected conformer policy.
Selection mismatch should be an error naming the first missing or different
identity.

The same selected pairs should define both the fit and the reported RMSD because
that matches the requested tool. A future API may deliberately separate
"alignment atoms" from "measurement atoms," as MDTraj permits by superposing
first and then evaluating without another fit
([MDTraj RMSD API](https://mdtraj.readthedocs.io/en/latest/api/generated/mdtraj.rmsd.html)),
but adding that flexibility now would complicate both terminology and output.

For the clustering workload, validate the shared chain/residue/atom identity
once against a reference structure, store coordinates in that canonical order,
and reuse them for all pairs. Re-parsing or rebuilding identity maps inside the
`O(S^2)` loop would obscure the solver benchmark and waste work.

## Community acceptance

Both approaches are mainstream:

- Biopython provides both an SVD superimposer and a QCP superimposer for protein
  and crystal structures
  ([SVD documentation](https://biopython.org/docs/latest/api/Bio.SVDSuperimposer.html),
  [QCP documentation](https://biopython.org/docs/latest/api/Bio.PDB.qcprot.html)).
- MDAnalysis uses QCP for minimum RMSD and optimal rotation
  ([official documentation](https://docs.mdanalysis.org/stable/documentation_pages/lib/qcprot.html)).
- MDTraj uses QCP for its optimized, parallel RMSD path
  ([official documentation](https://mdtraj.readthedocs.io/en/latest/api/generated/mdtraj.rmsd.html)).
- Rust molecular-analysis projects such as `groan_rs` and `molar` implement
  Kabsch/SVD, but through their own molecular system abstractions rather than a
  small coordinate-slice API
  ([groan_rs RMSD API](https://docs.rs/groan_rs/0.11.3/groan_rs/system/rmsd/index.html),
  [molar source](https://docs.rs/crate/molar/2.2.0/source/src/measure.rs)).

Thus "widely accepted" does not break the tie. QCP has stronger evidence for
high-throughput molecular RMSD; Kabsch has the broader general point-set and
linear-algebra footprint.

## Rust dependency survey

The survey used `cargo search` for `qcp`, `rmsd`, `kabsch`, and
`superposition`, followed by inspection of official crate metadata and source on
2026-08-28. Search results can change; the relevant primary indexes are
[crates.io `qcp`](https://crates.io/search?q=qcp),
[crates.io `rmsd`](https://crates.io/search?q=rmsd), and
[crates.io `kabsch`](https://crates.io/search?q=kabsch).

- `kabsch_umeyama` 0.1.2 is focused, but its public API uses const-generic array
  row counts, so a protein atom count must be known at compile time. It also
  pulls `nalgebra-lapack` while Arpeggia already has pure-Rust `nalgebra`, and it
  returns a general similarity transform rather than a protein-specific rigid
  fit. See its [official metadata](https://crates.io/crates/kabsch_umeyama/0.1.2)
  and [source](https://docs.rs/crate/kabsch_umeyama/0.1.2/source/src/lib.rs).
- `umeyama` 0.1.0 only accepts const-sized two-dimensional `f32` point arrays,
  so it is not applicable to protein structures
  ([official source](https://docs.rs/crate/umeyama/0.1.0/source/src/lib.rs)).
- `groan_rs` 0.11.3 supplies mass-weighted Kabsch through its full GROMACS-style
  `System` and group model and brings a large, overlapping dependency surface
  ([official crate metadata](https://crates.io/crates/groan_rs/0.11.3)).
- `molar` 2.2.0 similarly embeds Kabsch in a complete trajectory and molecular
  modeling library and is not a small numerical dependency
  ([official crate metadata](https://crates.io/crates/molar/2.2.0)).
- `cyanea-struct` 0.1.1 exposes a Kabsch result but introduces a separate protein
  representation and private linear-algebra implementation, duplicating both
  `pdbtbx` and `nalgebra` concerns already present here
  ([official crate metadata](https://crates.io/crates/cyanea-struct/0.1.1),
  [official source](https://docs.rs/crate/cyanea-struct/0.1.1/source/src/superposition.rs)).

The most economical Rust implementation is therefore a local coordinate-slice
kernel on top of the existing `nalgebra` dependency. If a future benchmark
reopens QCP, the original authors' ANSI C implementation is published under a
BSD license according to the 2010 paper, but any translation must retain
attribution and be checked against the repository's GPL-3.0 distribution
requirements
([Liu et al. 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2958452/)).

## Fit with the current repository

Arpeggia already depends on `nalgebra = 0.35`. The existing
`ResidueExt::center_and_normal` constructs a 3-by-N centered matrix and uses its
least singular vector as a fitted plane normal in
[`src/contacts/residues.rs`](../../src/contacts/residues.rs). Kabsch would also
use SVD, but on a different 3-by-3 cross-covariance matrix. There is no useful
shared "SVD algorithm" to extract beyond the dependency itself.

The implementation uses a small superposition module for coordinate centering,
proper-rotation validation, transform conventions, and RMSD. Moving `Plane`
there is unnecessary; a generic SVD abstraction shared only by plane fitting
and Kabsch would add indirection without removing meaningful duplication.

The transform convention must be documented and tested: which structure moves,
whether vectors are rows or columns, multiplication side, and whether translation
is applied before or after rotation. Existing libraries differ on these details;
for example, Biopython documents a right-multiplying rotation
([Biopython QCP API](https://biopython.org/docs/latest/api/Bio.PDB.qcprot.html)).

## Performance verification

Use representative structures and selection sizes (C-alpha, backbone, heavy
atoms, all atoms), and representative set sizes. Measure:

- preprocessing and identity validation once per input;
- total wall time and peak memory for the packed upper-triangle RMSD matrix;
- serial-pair and pair-parallel Kabsch throughput with identical cached
  coordinates;
- single-thread and normal Rayon settings;
- RMSD absolute differences between thread counts and maximum asymmetry;
- rotation determinant and transformed-coordinate residuals.

Do not add per-atom threading unless these end-to-end measurements show a
material gap on the supported workload. A single pair-parallel boundary is the
lean default.

## Resolved implementation scope

RMSD is uniformly weighted and requires at least three non-collinear selected
atoms in each of its two structures. The default selection is every
coordinate-observed protein residue in every chain, using C-alpha atoms.
Alternate conformers and models follow Arpeggia's existing preparation and
diagnostic policies; selected atom identities must match exactly. The rigid
transform and Kabsch solver remain implementation details.
