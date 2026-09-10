# Cluster exact-correspondence structures with Kabsch and k-medoids

Arpeggia uses a serial `f64` Kabsch kernel and a packed pairwise RMSD matrix
with k-medoids clustering. The public selection grammar, table schemas, input
formats, cache behavior, and measured performance are documented in the
[structure-clustering guide](https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md).

[ADR 0009](https://github.com/y1zhou/arpeggia/blob/master/docs/adr/0009-separate-sequence-correspondence-from-rmsd-evaluation.md) extends
two-structure RMSD with optional sequence correspondence and refinement; the
exact-correspondence collection behavior below remains unchanged.

## Correspondence and superposition

Selected atoms must correspond exactly after model and conformer selection.
Their identities include chain, author residue number, insertion code, residue
name, and atom name. A mismatch fails rather than silently intersecting atom
sets. Sequence/structural alignment and weighting remain deferred at the
correspondence boundary in [the RMSD selection module](https://github.com/y1zhou/arpeggia/blob/master/src/rmsd/selection.rs).

The Superposition Selection determines one proper rigid-body transform; the
RMSD Selection is evaluated with that transform without recentering or
refitting. Both residue selections independently default to all recognized
amino acids and share one atom preset. Superposition requires at least three
non-collinear atom pairs; evaluation requires at least one finite pair.
Reflection and physical scale fitting are prohibited. Weights are uniform,
and the residual sum is normalized by the evaluation atom count.

The transform stays private. `kabsch_rmsd` fits and evaluates the same arrays;
semantically equal parsed selections reuse this prepared fast path.

Kabsch uses the existing `nalgebra` SVD dependency. Coordinate normalization,
centering, scaled residual norms, and conditioning checks prevent overflow,
underflow, and misleading finite results for extreme coordinates. Numerical
normalization never fits a molecular scale factor. Fitted-versus-identity
fallbacks depend only on the Superposition Selection. An RMSD-only coordinate
that cannot be represented in that numerical frame produces a typed error;
Arpeggia does not maintain a second extreme-range scaling system. The detailed
numerical invariants and failure thresholds live beside the solver and tests.

Kabsch and QCP solve the same least-squares objective given identical paired
atoms, weights, and proper-rotation constraints. The solver is therefore an
implementation choice, not a public scientific-method option. Kabsch is short,
auditable, and needs no new dependency. QCP remains a possible optimization
only after an end-to-end benchmark shows a need.

Sources: [Kabsch 1976](https://doi.org/10.1107/S0567739476001873),
[Kabsch 1978](https://doi.org/10.1107/S0567739478001680),
[Theobald 2005](https://doi.org/10.1107/S0108767305015266), and
[Liu, Agrafiotis, and Theobald 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2958452/).
See the [superposition research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/structure-superposition.md) for
solver comparisons, QCP timing limitations, and dependency evidence. Plane fitting
uses SVD on a different matrix and does not justify a shared solver abstraction.

## Clustering objective and determinism

K-medoids consumes RMSD dissimilarities directly and selects observed structures
as representatives. It minimizes distances to medoids, not squared RMSDs or
coordinate distances in an invented common frame. Every observation receives a
cluster, matching the non-null cluster-ID contract. The existing GPL project
can use `kmedoids` 0.5.5 with default features disabled and an internal adapter
over Arpeggia's matrix.

A fixed count uses deterministic PAM BUILD initialization followed by FasterPAM.
A bounded automatic count uses DynMSC's medoid-silhouette objective. Automatic
selection is criterion-dependent, not parameter-free discovery of biological
states. Its maximum excludes `n`, because the all-singleton partition has a
trivially maximal silhouette. Numerically identical ensembles short-circuit to
one deterministic cluster. Fixed-count precedence and count domains are part
of the public guide.

Inputs are ordered by exact case-sensitive ID. Fixed-count equal-loss medoids
use the lowest canonical input index; a tie move is followed by reoptimization
because it can expose an improving swap. DynMSC retains its silhouette
objective instead of inheriting the different PAM-loss tie rule. Cluster labels
are normalized deterministically, and original Angstrom-valued distances remain
available for output.

Only the smallest uniform divisor needed to keep aggregate PAM losses finite
is applied. Positive uniform scaling preserves the objective; a distance range
that collapses relative to its maximum instead produces a calculation error.
The adapter avoids a square matrix copy and retains raw distances for ordinary
inputs. These numerical and convergence invariants are documented and tested
in [the clustering module](https://github.com/y1zhou/arpeggia/blob/master/src/clustering.rs).

Iteration exhaustion is a calculation failure. FasterPAM's cumulative swap
count sometimes requires a diagnostic pass to distinguish final-pass convergence
from exhaustion. DynMSC exposes only aggregate iterations across stages, so
complete stage-budget exhaustion can be detected but one isolated exhausted
stage cannot. Duplicate structures can also stop BUILD before the requested
maximum; the diagnostic uses the number actually produced.

Sources: [Kaufman and Rousseeuw 1987](https://wis.kuleuven.be/statdatascience/robust/papers/publications-1987/kaufmanrousseeuw-clusteringbymedoids-l1norm-1987.pdf),
[Schubert and Rousseeuw 2021](https://arxiv.org/abs/2008.05171),
[Lenssen and Schubert 2024](https://www.sciencedirect.com/science/article/pii/S0306437923001266),
and the [kmedoids adapter API](https://docs.rs/kmedoids/0.5.5/kmedoids/arrayadapter/index.html).

Average linkage remains a possible extension for a concrete hierarchy or
RMSD-cutoff requirement. Density, spectral, and affinity methods introduce
parameters or output semantics outside the fixed/automatic-count contract.
The [clustering research](https://github.com/y1zhou/arpeggia/blob/master/docs/research/structure-clustering.md) retains the
alternative-method comparison, additional crate survey, and supporting sources.

## Storage and execution boundaries

One packed `f64` triangle avoids a square copy. Each structure retains only the
union of fitting and evaluation coordinates, partitioned into superposition-only,
overlap, and RMSD-only blocks. Two shared boundaries describe the selections
without duplicated coordinates or indexed atom access in the hot loop. Raw
structures and non-reference identity tables are discarded after preparation.

`sysinfo` provides available-memory queries across supported platforms and
current-process cgroup limits on Linux. Only RAM and the needed process limits
are refreshed. Available RAM, rather than total or
merely free RAM, is the relevant estimate of reusable capacity. A matrix-only
preflight precedes parsing; a second check includes coordinate storage after
preparing the first structure. The guide records the estimates, exclusions,
80% threshold, unavailable-RAM warning, and explicit bypass. This is a heuristic,
not a guarantee that allocation will succeed. See the
[sysinfo memory API](https://docs.rs/sysinfo/0.39.6/sysinfo/struct.System.html#method.available_memory).

Preparation defines reference identities serially, then uses at most eight
parser workers to bound I/O and parser-memory pressure. Diagnostics return in
canonical input order. Pairwise workers own disjoint packed-matrix chunks;
each Kabsch atom loop and SVD remains serial. Each unordered pair is evaluated
once in canonical order, avoiding locks and schedule-dependent floating-point
reductions. Nested atom-level parallelism is deferred until measured need.

FasterPAM and DynMSC stay serial. Enabling parallel k-medoids would add
`ndarray` and randomization without providing parallel DynMSC. The
[local measurements](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/structure-clustering.md#local-structure-clustering-benchmark)
support keeping pairwise RMSD as the parallel boundary.

## Persistence and validation

Requested pairwise output is persisted before clustering to preserve useful
work if clustering fails. CLI cache reuse validates schema, complete pair
coverage, and exact IDs only; coordinate and selection provenance remain the
caller's responsibility. Malformed caches fail without overwrite, and new
caches use no-clobber creation.

Readers project required columns and reject wrong-size caches before complete
table materialization. The original lazy NDJSON reader was replaced by an eager
reader with a bounded row-count preflight after the
[v0.9.2 cleanup size and performance audit](https://github.com/y1zhou/arpeggia/blob/master/docs/research/v0.9.2-cleanup-audit.md).
CSV and Parquet use eager readers. XLSX remains excluded because it adds unrelated
reader/writer dependencies.

Equal-selection generalized results must match the prepared fast path bit for
bit in the same direction; reverse-direction and analytical checks use narrow
numerical tolerances. The independent-selection change retained its existing
coordinate payload for equal selections, with runtime gates of 5% for one
worker and 10% for eight, and a 10% peak-RSS gate. Overlapping selections must
save exactly `24n(f+r-u)` coordinate bytes. The guide retains the measurements;
regressions live with the implementation.

Exact correspondence, one shared atom preset, quadratic matrix storage,
heuristic memory protection, and caller-managed cache provenance remain the
principal limits. Broader polymer/ligand selection, alignment, weights, other
clustering methods, and public transforms require separate decisions.
