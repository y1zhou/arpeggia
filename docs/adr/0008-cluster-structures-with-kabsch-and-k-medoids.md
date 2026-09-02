# Cluster exact-correspondence structures with Kabsch and k-medoids

Arpeggia will add `rmsd`, `pairwise_rmsd`, and `cluster_structs` calculations
for ensembles whose structures contain the same protein chains and residue
composition. The first release requires exact selected-atom correspondence.
Sequence or structural alignment is deferred, but residue-selection code will
carry a `TODO:` identifying that future correspondence boundary.

## Inputs and selection

An ensemble comes from a non-recursive directory or a CSV, Parquet, or NDJSON
manifest. Directory extensions are matched case-insensitively for PDB
and mmCIF files; filename stems remain case-sensitive IDs. Manifests use `id`
and `path` by default, allow caller-selected column names, resolve relative
paths against the manifest directory, and ignore unrelated columns. XLSX is
excluded because it would add reader and writer dependencies.
Structure, manifest, and cache paths must resolve to regular files.

Before calculations begin, Arpeggia rejects fewer than two structures for
pairwise RMSD, fewer than three structures for clustering, empty or duplicate
IDs, duplicate canonical paths, missing paths, unsupported formats, and invalid
arguments. Inputs are sorted by exact case-sensitive ID for deterministic
pair ordering and cluster labels.

The Superposition Selection defaults to every coordinate-observed amino acid
recognized by Arpeggia in every chain and C-alpha atoms. The RMSD Selection
independently defaults to all recognized residues. Callers use the same residue
string for both selections when they want to superpose and evaluate a subset.
Both selections use one shared atom preset: `ca`, `backbone`, `heavy`, or `all`;
backbone means `N`, `CA`, `C`, `O`, and `OXT`. Heavy mode excludes hydrogen and
deuterium, including digit-leading atom names when element metadata is absent.
Heavy and all-atom modes include selected ACE/NH2 terminal caps but exclude
solvent, ions, and ligands.
Both selections use comma-separated author-numbered chain/residue clauses such as
`A`, `A:1-100`, `A:110-120`, `B`, `C:1`, `C:3`, and `C:9-20`. Insertion codes
and negative author residue numbers are valid, so `A:-5--1` selects residues
-5 through -1. Overlaps are deduplicated; malformed clauses, unknown chains,
and empty selections fail. Both grammars are validated before structure I/O;
later chain, correspondence, and nonempty-selection errors identify whether the
Superposition Selection or RMSD Selection failed.

Model zero selects the first coordinate model. An explicitly selected model is
applied uniformly to the ensemble. Automatic first-model or conformer selection
uses Arpeggia's existing structured diagnostics. Selected atoms are paired by
chain, author residue number, insertion code, residue name, and atom name after
conformer selection. Any mismatch fails rather than silently intersecting atom
sets.

## Superposition and RMSD

The numerical kernel accepts two equal-shaped `(n, 3)` `f64` coordinate arrays,
uses uniform atom weights, and returns an `f64` RMSD. It implements Kabsch with
the existing `nalgebra` SVD, forbids reflection, and does not estimate scale.
Each array is first expressed as anchor-relative displacements and normalized
by its largest displacement. If an anchor subtraction would overflow, raw
coordinates are scaled before that subtraction instead. The kernel caches the
normalized centroid, mean-centers coordinates while forming covariance, and
forms residuals before centering them so small orthogonal differences are not
rounded out of large coordinates. It expresses both arrays against one shared
physical magnitude and restores the RMSD to Angstroms only after a scaled norm.
The unrotated residual is also a numerical upper bound when SVD noise in a
near-degenerate covariance makes the fitted rotation worse. Its first-order
rotational component is checked separately so an unresolved small rotation is
not mistaken for physical deformation. This
translation-invariant normalization prevents centering, squaring, and solver
overflow or underflow without fitting a physical scale factor. The
non-collinearity precheck likewise normalizes coordinates before forming
displacements and cross products, so finite very large and very small
structures do not overflow or underflow that validation.
If either fitted residual uncertainty or an unresolved first-order rotation is
solver-scale in normalized units but restoring physical units would amplify it
above `1e-6` Angstrom, the kernel returns a typed calculation error instead of
a misleading finite RMSD.
At least three non-collinear Superposition Selection atoms are required. The
RMSD Selection requires at least one corresponding finite atom pair and need
not be non-collinear. Kabsch determines one proper rigid-body transform solely
from the Superposition Selection; that transform is then applied to the RMSD
Selection without recentering or refitting it. The rigid transform and solver
remain private. The public coordinate-array `kabsch_rmsd` convenience function
continues to use the same arrays for superposition and RMSD evaluation. Any
fitted-versus-identity numerical fallback is chosen solely by the
Superposition Selection objective. RMSD uses uniform atom weights and divides
the squared residual sum by the RMSD Selection atom count. The complete union
is expressed in the Superposition Selection's numerical frame. If a
pathological RMSD-only coordinate cannot be represented finitely in that frame,
the calculation returns a typed error rather than introduce a second scaling
system for physically meaningless coordinate ranges.

When the normalized parsed selectors are semantically equal, the production
path reuses the existing prepared Kabsch RMSD calculation without a redundant
evaluation pass. Tests and local benchmarks must also exercise the generalized
transform-plus-evaluation path with equal selections and compare its result
against the existing calculation; the fast path must not conceal a numerical
disagreement.

Kabsch is selected over QCP because it is short, auditable, broadly accepted,
and needs no new superposition dependency. QCP solves the same least-squares
objective and may be faster for some workloads, but maintaining a local second
solver is not justified without measured end-to-end need. The existing plane
fit also uses SVD, but its matrix and purpose differ; no generic SVD abstraction
will be introduced.

The scalar Rust/Python `rmsd` interface and `rmsd` CLI command take exactly two
structures, accept `superpose_residues` and `rmsd_residues`, and remain serial.
The CLI arguments are `--superpose-residues`/`-s` and
`--rmsd-residues`/`-r`; both are ordinary strings that independently default to
the existing empty-string selection of all recognized residues. The CLI prints
the scalar value in Angstroms to stdout.
Python and Rust also expose `pairwise_rmsd` for an ensemble. Python
returns a Polars DataFrame in long form with exactly `id_1`, `id_2`, and `rmsd`;
it contains each unordered pair once and excludes the diagonal. Scalar and
pairwise outputs report only the RMSD Selection result, not a second
Superposition Selection RMSD.

## Clustering

`ClusteringMethod` initially contains only `KMedoids`. The implementation uses
`kmedoids` 0.5.5 with default features disabled. A fixed `num_clusters` uses
deterministic PAM BUILD initialization followed by FasterPAM. A supplied
`max_clusters` uses DynMSC to select from 2 through that maximum. One of the
two bounds is required; when both are present, `num_clusters` wins and a
structured warning reports that `max_clusters` was ignored. Fixed counts allow
`1..=n`; automatic maxima allow `2..=(n-1)`. Allowing `max_clusters=n` would make
the all-singleton partition achieve the maximum medoid-silhouette score by
construction. An all-zero matrix under automatic
selection short-circuits to one deterministic cluster; RMSDs at or below
`1e-12` Angstrom are numerical zero for this check. The default iteration limit
is 100. FasterPAM limit exhaustion is a Calculation Failure. When the run uses
the full iteration count and made at least one swap, one separate diagnostic
pass accepts the result if no further improving swap exists and otherwise
reports exhaustion. The upstream API returns only cumulative swaps, so this
diagnostic can be redundant when the final allowed pass itself made no swap; a
run with no swaps at all is already converged. The upstream
DynMSC API exposes only aggregate iterations across its internal cluster-count
stages, so Arpeggia can detect complete stage-budget exhaustion but cannot
distinguish one isolated stage reaching its limit. The aggregate ceiling uses
the number of medoids PAM BUILD actually produced because duplicate structures
can stop BUILD before the requested maximum.

For fixed `k > 1`, equal-loss within-cluster medoid ties move to the lowest
canonical input index. Because a plateau move can expose a subsequent
improving swap, FasterPAM resumes within the remaining iteration budget after
medoid identities and canonical memberships stabilize; exhaustion is diagnosed
as above. PAM BUILD
already resolves `k=1` ties by first input order, so that case avoids the
redundant quadratic tie pass. DynMSC retains its own deterministic
medoid-silhouette result rather than applying the different PAM-loss tie rule.

Before invoking PAM, matrix lookups use the smallest uniform divisor needed to
keep worst-case aggregate losses below the finite `f64` range. Uniform positive
scaling cannot change the PAM objective or partition, and using the smallest
divisor preserves mixed-scale distances. If a positive distance vanishes when
expressed relative to the matrix maximum, clustering returns a typed calculation
error rather than optimize a collapsed distance matrix.
Original Angstrom-valued RMSDs remain stored for output and `rmsd_to_medoid`.
Ordinary matrices needing no divisor return raw distances from the adapter;
extreme matrices multiply by one cached reciprocal in the PAM hot path.

The cluster table contains `id: String`, `cluster_id: UInt32`,
`medoid_id: String`, and `rmsd_to_medoid: Float64`. Cluster IDs are contiguous
from zero and normalized deterministically by canonical input order. A medoid
is an observed structure, not an averaged centroid.

Python `cluster_structs` requires exactly one of a structure input or an
explicit `pairwise_rmsd` Polars DataFrame. A supplied table must contain at
least three IDs, every unordered pair exactly once, no diagonal, finite
non-negative numeric RMSDs, and string IDs; numeric RMSDs are converted to
`Float64`, and unrelated columns are ignored. The CLI always starts from a
non-recursive structure directory and does not accept a manifest or arbitrary
matrix input.

## Memory and parallel execution

The internal matrix stores one `f64` triangle without diagonal, using
`4n(n-1)` bytes. Before parsing structures, Arpeggia compares that estimate
with 80% of effective available RAM. Effective availability is host available
RAM capped by current-process cgroup availability when reported by `sysinfo`.

If the matrix passes, Arpeggia prepares the first structure and calculates a
second estimate for the matrix plus `24nu` bytes of selected `(x, y, z) f64
coordinates, where `u` is the union atom count across both selections. This
check occurs before loading the remaining structures. The
estimate excludes parser transients, identity keys, allocator overhead, output
DataFrames, and clustering scratch space and is documented as a heuristic, not
an allocation guarantee. If availability cannot be queried, an estimate above
8 GiB produces a warning instead of a hard limit. `bypass_mem_check` and
`--bypass-mem-check` skip both checks, warnings, and failures.

Only coordinates in the union of the Superposition Selection and RMSD Selection
are retained. Within each structure they are ordered as superposition-only,
overlap, then RMSD-only atoms. The two selections are therefore contiguous,
possibly overlapping slices described by two shared boundaries, without
coordinate duplication, masks, or indexed access in the pairwise hot loop.
One atom traversal constructs identities and classifies union members; the
temporary classification is discarded after partitioning. Each full parsed
structure is dropped after preparation. The first structure is prepared
serially to define
the reference identities and atom counts. Each later atom-identity table is
discarded immediately after Exact Atom Correspondence is validated for both
selections. Future sequence or structural correspondence can produce a compact
mapping at this boundary without retaining all raw identity tables. Remaining
structures use
`min(num_threads, 8)` parser workers, with automatic thread selection capped at
eight. Successful results and diagnostics are restored to canonical input
order; an error cooperatively stops collection of later parsing work.

Pairwise matrix construction then uses up to the smallest of the requested
Rayon worker count, available processors, and the number of structure pairs.
Workers own disjoint contiguous packed-matrix chunks and run the Kabsch atom loop and
3-by-3 SVD serially for each pair. This avoids locks, nested parallelism, and
schedule-dependent floating-point reductions. Each unordered pair is evaluated
once in canonical ID order; reverse-direction equivalence is tested within a
strict numerical tolerance rather than doubling work and averaging two fits.
The current `nalgebra` and
`matrixmultiply` feature graph is single-threaded for this operation, so no
manifest change is needed. FasterPAM and DynMSC remain serial; parallel
k-medoids adds `ndarray`, randomization, and no parallel DynMSC implementation,
and will be reconsidered only if later benchmarks identify clustering as a
serious bottleneck.

`num_threads` controls parser preparation and pairwise matrix construction in
`pairwise_rmsd` and structure-backed `cluster_structs`. It is omitted from the
scalar `rmsd` operation and has no effect when Python clustering receives an
already calculated DataFrame.

## Output ordering and cache reuse

Cluster and optional pairwise tables support CSV, Parquet, and NDJSON.
Cached Parquet row counts are checked from file metadata before columns are
decoded; cached CSV and NDJSON reads request at most one row beyond the expected
count. NDJSON projection and row limiting happen in one lazy Polars collection.
When pairwise output is requested, it must be written successfully before
k-medoids begins, preserving useful work if clustering subsequently fails.

On a later CLI run with pairwise output requested, Arpeggia checks only the
exact requested output path. A present table is reused when its schema,
complete unordered-pair coverage, and exact ID set validate. Coordinate files,
superposition and RMSD selection settings, model choice, conformers, and
Arpeggia version are
deliberately not checked. Documentation warns that cache validity is therefore
the caller's responsibility, and a debug log names the reused file. Removing
the file forces recalculation. A malformed, incomplete, or ID-mismatched file
causes an error with removal instructions; it is never silently overwritten.
New caches use no-clobber creation so a concurrently created path is preserved.

## Consequences

This design favors a small deterministic scientific core and an existing
specialized clustering crate. Its main limits are exact atom correspondence,
one shared atom preset across the two residue selections, quadratic time and
matrix storage, heuristic rather than guaranteed memory
protection, the recognized-amino-acid selection domain, and caller-managed CLI
cache provenance. Future exact-correspondence all-atom selection may retain
arbitrary polymers, ligands, and modified residues when both inputs match.
Sequence or structural correspondence followed by weighted RMSD, alternate
clustering algorithms, per-atom RMSD parallelism, parallel k-medoids, and public
rigid transforms require separate evidence and decisions.

Research and source comparisons are recorded in
[`structure-superposition.md`](../research/structure-superposition.md) and
[`structure-clustering.md`](../research/structure-clustering.md).

Completion requires regression tests for equal, disjoint, partially overlapping,
and one-atom RMSD Selections; a synthetic rigid-domain movement superposed on one
chain and evaluated on others; selection-specific failures; independent defaults;
and forward/reverse equivalence. Local benchmarks compare the existing and new
implementations for equal, heavily overlapping, and disjoint selections,
recording exact retained coordinate payload and observed peak RSS. The retained
payload compares `24nu` with the duplicated-array counterfactual `24n(f+r)`.
For equal selections, median runtime must remain within 5% of the existing
single-worker path and within 10% with eight workers; exceeding either threshold
requires investigation. These are runtime and memory gates, not tolerances for
calculated RMSD values. The generalized equal-selection calculation must match
the existing prepared Kabsch result bit for bit in the same direction; only
reverse-direction symmetry and independently derived analytical expectations
use numerical tolerances. Equal selections must retain exactly the existing
coordinate payload and may not increase median peak RSS by more than 10%.
Overlapping selections must demonstrate the exact payload saving
`24n(f+r-u)` predicted by union storage.
