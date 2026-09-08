# Protein structure superposition and RMSD

Research date: 2026-08-28

This note retains solver comparisons, external performance evidence, and the
Rust dependency survey. The accepted behavior is maintained in
[ADR 0008](../adr/0008-cluster-structures-with-kabsch-and-k-medoids.md):

- [Correspondence and superposition](../adr/0008-cluster-structures-with-kabsch-and-k-medoids.md#correspondence-and-superposition): exact atom pairing, independent fit/evaluation selections, numerical safeguards, and the choice of Kabsch.
- [Storage and execution boundaries](../adr/0008-cluster-structures-with-kabsch-and-k-medoids.md#storage-and-execution-boundaries): coordinate preparation, matrix storage, and pair-level parallelism.
- [Persistence and validation](../adr/0008-cluster-structures-with-kabsch-and-k-medoids.md#persistence-and-validation): numerical equivalence and performance gates.

Current selection syntax and examples belong in the
[usage guide](../benchmarks/structure-clustering.md#rmsd), and the numerical implementation
and regressions live in [`src/rmsd.rs`](../../src/rmsd.rs).

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

## QCP numerical evidence

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
