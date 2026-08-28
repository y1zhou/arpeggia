# Clustering a pairwise protein-structure RMSD matrix

Date: 2026-08-28

## Scope and recommendation

This note evaluates clustering after Arpeggia has already computed one
least-squares-superposed RMSD for every pair of structures. It assumes that
each structure contributes the same ordered atom set. Atom matching,
superposition, and the later sequence-alignment extension are separate design
questions.

For the first implementation, use **k-medoids** as the primary method:

- for a requested cluster count, initialize deterministically with PAM BUILD
  and optimize with FasterPAM;
- for an automatically selected count, use DynMSC over an explicit
  `k_min..=k_max` range;
- report an actual input structure as each cluster's **medoid**, plus every
  member's RMSD to that medoid.

The fit is unusually good for this problem. PAM is defined on pairwise
dissimilarities rather than coordinate vectors, and its center is the observed
object with the smallest dissimilarity to the other members. Thus it consumes
the RMSD matrix directly and its representative is a real, writable protein
structure rather than an averaged coordinate model. FasterPAM retains PAM's
result while improving its SWAP phase by a factor of order `k`; eager swaps
trade exact swap order for comparable quality and more speed
([Schubert and Rousseeuw, 2021](https://arxiv.org/abs/2008.05171)). This center
definition also matches the one used by GROMACS for structure clustering: its
cluster center is the structure with the smallest average RMSD to the other
members
([GROMACS `gmx cluster`](https://manual.gromacs.org/current/onlinehelp/gmx-cluster.html)).

Average-linkage agglomerative clustering remains a future option if a
hierarchical view or RMSD cutoff becomes a concrete release requirement. It is
a defensible method with a small Rust dependency, but it adds method and cutoff
choices without improving the accepted fixed/automatic-k workflow.

Do not expose DBSCAN, HDBSCAN, spectral clustering, or affinity propagation in
the first release. They have legitimate later uses, but their output or
parameters conflict with the proposed simple API, and the available Rust
implementations are a poorer match for an existing dense dissimilarity matrix.

## Why k-medoids is the best default

The original PAM objective chooses `k` representative observations and
minimizes the total dissimilarity from every observation to its nearest
representative
([Kaufman and Rousseeuw, 1987](https://wis.kuleuven.be/statdatascience/robust/papers/publications-1987/kaufmanrousseeuw-clusteringbymedoids-l1norm-1987.pdf)).
For Arpeggia, use RMSD itself as the dissimilarity. Squaring RMSD would define a
different, more outlier-sensitive objective and should not happen implicitly.

Advantages for this use case:

- No Euclidean feature embedding or distance-to-similarity transform is
  required. The implementation supports arbitrary dissimilarity matrices
  ([`rust-kmedoids` documentation](https://github.com/kno10/rust-kmedoids)).
- Every structure is assigned to exactly one cluster. This preserves the
  requested non-null `UInt32` `cluster_id` naturally.
- A medoid is an observed structure, so `medoid_id` and `rmsd_to_medoid` have
  direct physical and operational meanings.
- The matrix is already the expensive and unavoidable part of the requested
  workflow. K-medoids can read the same condensed matrix without inventing a
  coordinate embedding.
- A maintained, domain-authored Rust implementation exists. `kmedoids` 0.5.5
  implements PAM, FasterPAM, parallel FasterPAM, silhouette variants, and
  DynMSC, and accepts either an `ndarray::Array2` or its own serialized
  `LowerTriangle` adapter
  ([crate API](https://docs.rs/kmedoids/latest/kmedoids/),
  [array adapters](https://docs.rs/kmedoids/latest/kmedoids/arrayadapter/index.html)).
  Its license is GPL-3.0-or-later, compatible with this GPL project
  ([crate manifest](https://docs.rs/crate/kmedoids/latest/source/Cargo.toml.orig)).

### Fixed `k`

Use `pam_build(matrix, k)` to obtain deterministic initial medoids, then pass
those medoids to `fasterpam`. The crate documents BUILD as an initializer and
FasterPAM as the higher-quality optimization step
([`pam_build`](https://docs.rs/kmedoids/latest/kmedoids/fn.pam_build.html),
[`fasterpam`](https://docs.rs/kmedoids/latest/kmedoids/fn.fasterpam.html)).
With structures sorted into a canonical input order and an explicit tie rule,
this avoids seed-dependent output. Random restarts may find another local
optimum, but they would make a stable scientific CLI harder to reproduce and
are not justified before benchmarks show a material quality gain.

The crate says its parallel implementation is typically faster above roughly
5,000 observations. Enabling its `parallel` feature also enables `rayon`,
`rand`, and `ndarray`; Arpeggia already has Rayon, while it does not otherwise
depend directly on `ndarray`
([`rust-kmedoids` README](https://github.com/kno10/rust-kmedoids),
[crate manifest](https://docs.rs/crate/kmedoids/latest/source/Cargo.toml.orig)).
The lean initial dependency is therefore:

```toml
kmedoids = { version = "0.5.5", default-features = false }
```

Parallelize the much more expensive all-pairs RMSD calculation with Arpeggia's
existing Rayon pool. Revisit parallel FasterPAM only after an end-to-end
benchmark above 5,000 structures.

### Automatic `k`

DynMSC begins at a supplied maximum number of medoids, directly optimizes the
average medoid silhouette, removes one medoid at a time, and retains the best
solution down to a supplied minimum
([DynMSC API](https://docs.rs/kmedoids/latest/kmedoids/fn.dynmsc.html)). The
underlying publication presents it specifically as medoid-silhouette
clustering with automatic cluster-number selection
([Lenssen and Schubert, 2024](https://www.sciencedirect.com/science/article/pii/S0306437923001266)).

This is "automatic within declared bounds," not a parameter-free discovery of
ground truth. Internal cluster validation is inherently criterion-dependent;
the DynMSC paper explicitly frames cluster evaluation as subjective and
data-dependent. The accepted design fixes `k_min` at 2 and requires `k_max`,
then returns the selected `k` in diagnostics. The Python package maintained by
the same authors warns that an excessively high `k_max` can favor many
singleton clusters and recommends a maximum only two to three times the
expected count
([official Python wrapper README](https://github.com/kno10/python-kmedoids)).

Recommended edge behavior:

- reject clustering inputs with fewer than three structures; a singleton has
  no clustering interpretation, and a pair does not justify this interface;
- require at least two structures for an ensemble pairwise-RMSD calculation;
- all pairwise RMSDs zero: one cluster, with a deterministic medoid;
- otherwise default `k_min >= 2`, because the usual silhouette comparison
  needs another cluster;
- reject `k_max > n` and non-finite or negative dissimilarities before calling
  the crate, whose low-level functions document panics for invalid `k`.

## Algorithm comparison

| Method | Uses RMSD matrix directly? | Cluster count | Representative | Scientific/API fit | Runtime and memory after RMSDs | Rust status |
|---|---|---|---|---|---|---|
| PAM / FasterPAM | Yes; arbitrary dissimilarity | Fixed `k` | Medoid, an input structure | Best default; all observations assigned and objective is transparent | Dense/condensed matrix remains `O(n^2)`; FasterPAM gives an order-`k` SWAP speedup | `kmedoids` 0.5.5 is a direct fit |
| DynMSC | Yes | Chooses within `k_min..=k_max` | Medoid | Best available automatic option for the same model; bounds and criterion must be disclosed | Reuses the matrix and optimizes several decreasing `k` values | Included in `kmedoids` 0.5.5 |
| Agglomerative, average linkage | Yes, condensed | Cut dendrogram to fixed `k` or at an RMSD threshold | None intrinsically; compute a post-hoc medoid | Good optional hierarchy; all observations assigned | Average linkage is `O(n^2)` time and memory | `kodama` 0.3.0 consumes condensed dissimilarities directly |
| DBSCAN | Conceptually yes | Emerges from `eps` and `min_samples` | None | Useful for dense states plus noise, but `eps` is effectively a scientific RMSD cutoff and some points are noise | A dense precomputed matrix gives `O(n^2)` storage/search rather than the spatial-index advantage | Linfa works on feature rows; no compelling direct precomputed adapter found |
| HDBSCAN | Conceptually yes | Automatic from density hierarchy and minimum cluster size | None; returns outlier scores | Better than DBSCAN for variable density, but still permits noise and adds density semantics | Dense fallback is at least quadratic; spatial trees cannot exploit a precomputed opaque matrix automatically | `petal-clustering` 0.13.0 requires point rows and constructs a `BallTree`; not a direct fit |
| Spectral | Requires an RMSD-to-affinity kernel | Normally fixed `k`; eigengap is another heuristic | None; post-hoc medoid | Can find non-convex graph partitions, but adds an arbitrary kernel scale and loses the simple RMSD objective | Dense eigendecomposition is approximately `O(n^3)` and the affinity is `O(n^2)` | Available crates target feature/graph input, not this lean pipeline |
| Affinity propagation | Yes after choosing a similarity such as `-RMSD` or `-RMSD^2` | Emerges from a preference parameter | Exemplar, an input structure | Attractive automatic exemplars, and used by MDAnalysis, but preference controls `k` indirectly and convergence is not guaranteed | `O(n^2)` time and memory per documented Rust implementation, with several dense message arrays | `affinityprop` 0.2.0 accepts precalculated similarity, but adds `ndarray` and more state |

### Agglomerative clustering details

Average linkage defines inter-cluster distance as the mean of every
cross-cluster pairwise distance. Complete linkage uses the maximum, and single
linkage the minimum
([SciPy linkage definitions](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)).
From those definitions, average linkage is the most neutral first choice for
an RMSD ensemble: it does not let one nearest pair chain whole conformational
states together as single linkage can, and it does not let the single most
distant pair dictate every merge as complete linkage does. That suitability
statement is an inference from the linkage definitions, not a universal claim
that average linkage recovers a biological ground truth.

Do not offer Ward, centroid, or median linkage over the superposition RMSD
matrix. Their Lance-Williams updates are correctly defined only for Euclidean
pairwise distances in a common feature space; SciPy places that verification
responsibility on callers supplying a precomputed matrix
([SciPy linkage notes](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)).
Pairwise optimally superposed structures do not share one common Cartesian
frame, so Arpeggia should not silently claim that prerequisite.

`kodama` performs hierarchical clustering directly on a condensed
dissimilarity vector. Its average-linkage implementation is `O(n^2)`, and its
runtime is documented as comparable to `fastcluster`
([`kodama` documentation](https://docs.rs/kodama/latest/kodama/)). A thin
Arpeggia layer would still need to cut its dendrogram and calculate each
cluster's medoid. Caveat: `kodama` explicitly says tie behavior is unspecified
and that comparison with SciPy/fastcluster was hand-checked rather than
automated. If exposed, Arpeggia must sort inputs canonically, pin the crate,
test tied matrices, and document that equally valid dendrograms may differ at
ties.

### Density-based methods

DBSCAN was designed to discover arbitrarily shaped dense regions and noise
using a neighborhood radius and minimum density
([Ester et al., 1996](https://file.biolab.si/papers/1996-DBSCAN-KDD.pdf)).
That is useful if rare, isolated conformers should be marked as outliers rather
than forced into a cluster. It is not parameter free, however: the RMSD radius
has a strong domain meaning, and its noise output conflicts with a required
non-null `cluster_id: UInt32`. Before adding DBSCAN, decide whether noise has a
nullable cluster ID, a reserved ID, or singleton clusters.

HDBSCAN replaces one global density cut with a hierarchy and selects stable
clusters, returning outliers and outlier scores
([Campello et al. reference and API in `petal-clustering`](https://docs.rs/petal-clustering/latest/petal_clustering/struct.HDbscan.html)).
The available Rust implementation is not wired for a precomputed matrix: its
source accepts an `ndarray` of feature rows, copies it to standard layout,
builds a `BallTree`, and queries the supplied metric
([`petal-clustering` source](https://docs.rs/petal-clustering/latest/src/petal_clustering/hdbscan.rs.html)).
Encoding row indices as fake feature vectors merely to recover distances from
an external matrix would defeat that API and its spatial tree. Defer until a
tested direct-distance implementation exists or there is a concrete outlier
requirement worth a small internal adapter.

### Spectral clustering and affinity propagation

The canonical spectral algorithm eigendecomposes a normalized affinity matrix
and then applies k-means to the selected eigenvectors
([Ng, Jordan, and Weiss, 2001](https://www.ee.columbia.edu/~dpwe/papers/NgJW01-specclus.pdf)).
RMSDs must first become affinities, typically through a scale-bearing kernel;
that extra scale can change the clusters. Dense eigendecomposition and a
second clustering algorithm are disproportionate when the direct
dissimilarity methods already match the desired output.

Affinity propagation accepts pairwise similarities, considers all points as
candidate exemplars, and uses self-preference values rather than a prescribed
cluster count
([Frey and Dueck, 2007](https://people.csail.mit.edu/kjhsiao/Frey2007.pdf)).
It is established in protein-ensemble software: MDAnalysis accepts a
triangular pairwise similarity matrix and exposes affinity propagation and
DBSCAN for conformational ensembles
([MDAnalysis/mdaencore clustering](https://www.mdanalysis.org/mdaencore/encore/clustering.html)).
But the original paper also documents oscillation/non-convergence under
degenerate similarities, mitigated by tiny noise or extra damping. Adding
noise undermines deterministic behavior, while preference and damping broaden
the public API. The Rust `affinityprop` crate documents precalculated
similarities, Rayon parallelism, and quadratic time/memory
([crate documentation](https://docs.rs/affinityprop)). It is a reasonable
future experiment, not the lean first choice.

## Matrix representation and scale

All exact methods considered here inherit the all-pairs bottleneck. Store one
triangle, excluding the zero diagonal, rather than an `n x n` square. The
accepted design uses `f64`, so the condensed matrix uses `4 * n * (n - 1)`
bytes: about 400 MB for 10,000 structures and 10 GB for 50,000. A square matrix
doubles those figures before algorithm scratch space.

Implementation implications:

- Calculate each pair once and store it in canonical condensed order.
- Do not build a Polars pairwise-results DataFrame unless the user requested
  that output. A long table has `n(n-1)/2` rows and substantial string/column
  overhead beyond the numeric matrix.
- If pairwise output is requested, derive or stream `id_1`, `id_2`, and `rmsd`
  from the matrix; do not recalculate RMSDs.
- Write requested pairwise output successfully before starting k-medoids. On a
  later CLI run, reuse an existing table at the exact requested output path
  when its validated ID set matches the current inputs. This intentionally does
  not verify coordinates or selection settings; document that limitation and
  emit a debug message when reuse occurs. Removing the file forces
  recalculation.
- Use that same three-column long form when a caller supplies a previously
  calculated matrix to the Python clustering API. A downstream pivot can
  construct a square representation when needed.
- Keep accepted tabular formats to CSV, Parquet, JSON, and NDJSON. Polars'
  existing `json` feature provides both JSON variants, so NDJSON requires no
  additional dependency or feature. XLSX would require another dependency and
  is deliberately outside the feature.
- Keep the clustering dependency behind an internal adapter so the RMSD
  engine, Python API, and CLI do not expose crate-specific types.
- Validate finite, symmetric, non-negative results and an exact zero diagonal
  at the matrix boundary.
- Benchmark end-to-end time by both structure count `n` and selected atom
  count `a`. RMSD work is `O(n^2 a)` and is likely to dominate clustering for
  realistic proteins.

The k-medoids crate's `LowerTriangle` adapter can avoid `ndarray`; `kodama`
expects a condensed vector and mutates it while constructing the dendrogram.
If both methods are exposed, define one Arpeggia-owned condensed matrix and
adapt at the method boundary. Do not retain simultaneous square, lower-triangle,
and output-table copies.

## Output semantics

Do not call the representative a centroid. A centroid is normally an average
coordinate vector and may not correspond to an input structure; after
independent pairwise fits there is not even one shared coordinate frame in
which that average is already defined. Use this stable core schema:

| Column | Type | Meaning |
|---|---|---|
| `id` | `Utf8` | Exact input ID |
| `cluster_id` | `UInt32` | Stable contiguous label `0..k-1` |
| `medoid_id` | `Utf8` | ID of the observed representative structure |
| `rmsd_to_medoid` | `Float64` | Selected-atom RMSD already present in the pairwise matrix |

For deterministic labels, sort directory and table inputs into a canonical
order before RMSD calculation. Then sort final clusters by medoid input index
(or, after deciding ID normalization, by `medoid_id`) and remap them to
contiguous `UInt32` IDs. In a distance tie, choose the lowest canonical input
index as medoid. The method's objective value, selected `k`, bounds, iteration
limit, and convergence status belong in structured diagnostics or run
metadata, not duplicated in every output row.

## Verification and benchmarks required before stabilizing the API

1. Synthetic rigid transformations must have RMSD approximately zero and
   remain in the same cluster.
2. Synthetic conformational states with known within/between RMSD separation
   should recover the expected partition for fixed-k FasterPAM, automatic
   DynMSC, and any exposed hierarchical method.
3. Permuting input rows must preserve the partition and medoid IDs after label
   normalization, except where equal-distance ties admit multiple documented
   optima.
4. Duplicate and all-identical structures, `n = 1`, `k = 1`, `k = n`, invalid
   bounds, non-finite coordinates/RMSDs, and repeated IDs require explicit
   tests.
5. Compare the fixed-k objective and medoids against the crate's own PAM on
   small matrices. If hierarchical clustering is exposed, compare partitions
   against SciPy average linkage for matrices with and without ties; `kodama`
   itself identifies tie behavior and automated cross-implementation testing
   as caveats.
6. Benchmark release builds at a grid of structure counts and atom counts,
   reporting pairwise-RMSD time, clustering time, peak memory, and optional
   pairwise-output cost separately. Include at least one case above the
   `rust-kmedoids` project's approximate 5,000-item threshold before enabling
   its parallel feature.

## Follow-up: crate availability and RAM limits

Date checked: 2026-08-28

### Off-the-shelf clustering crates

Average-linkage clustering does have a direct Rust implementation. `kodama`
0.3.0 accepts a mutable condensed pairwise dissimilarity vector and
`Method::Average`; its documentation gives average linkage `O(n^2)` runtime
and returns the complete dendrogram
([`kodama` API](https://docs.rs/kodama/0.3.0/kodama/)). The repository is not
archived, but development is quiet: its last push was 2025-04-09 when checked
through the GitHub repository API. `linfa-hierarchical` is maintained as part
of Linfa, but it delegates its agglomeration to `kodama` and presents a
similarity-kernel/`ndarray` interface instead of improving the fit for an
existing condensed RMSD matrix
([Linfa hierarchical documentation](https://docs.rs/linfa-hierarchical/0.8.1/linfa_hierarchical/)).

No maintained crate directly implementing the Daura/GROMOS conformation
algorithm over a precomputed dissimilarity matrix was found in the crates.io
index under either
[`gromos`](https://crates.io/search?q=gromos) or
[`daura`](https://crates.io/search?q=daura). This is distinct from the
algorithm being unavailable generally: GROMOS++ documents its `cluster`
program as consuming an RMSD matrix, repeatedly choosing the structure with
the most neighbors below a cutoff, and removing that cluster
([GROMOS++ manual](https://www.gromos.net/gromos11_pdf_manuals/vol5.pdf)).
Therefore, under the design rule that every exposed method must have a suitable
off-the-shelf Rust implementation, the first release should expose
**k-medoids only**.

### Querying available RAM

Rust's standard library does not expose host RAM statistics;
[`std::alloc::System`](https://doc.rust-lang.org/std/alloc/struct.System.html)
is an allocator, not a memory-capacity query. Of the external choices:

- `procfs::Meminfo` exposes Linux `MemAvailable`, defined as an estimate of
  memory available to start applications without swapping, but the crate is
  explicitly Linux-only and cgroup handling would remain Arpeggia's job
  ([`procfs` crate](https://docs.rs/procfs/0.18.0/procfs/),
  [`Meminfo`](https://docs.rs/procfs/0.18.0/procfs/struct.Meminfo.html)).
- `systemstat` supports Linux, Windows, and macOS, but its common `Memory`
  result exposes `total` and `free`, not a documented cross-platform
  `available` value. It also supplies many unrelated system-statistics APIs
  ([`systemstat::Platform`](https://docs.rs/systemstat/0.2.7/systemstat/platform/common/trait.Platform.html)).
- `sysinfo` supports all three target systems and directly exposes
  `available_memory()` in bytes after a RAM refresh. It distinguishes free
  (unallocated) memory from available (available for reuse) memory and provides
  Linux-only cgroup limits for both the system and a process
  ([`System` memory API](https://docs.rs/sysinfo/0.39.6/sysinfo/struct.System.html#method.available_memory),
  [`Process::cgroup_limits`](https://docs.rs/sysinfo/0.39.6/sysinfo/struct.Process.html#method.cgroup_limits)).

Use `sysinfo` and refresh RAM only; do not call `System::new_all()`. Arpeggia's
current Polars dependency already resolves `sysinfo` 0.39.6, so a compatible
direct dependency does not add another compiled crate version. Define
**effective available RAM** at startup as host `available_memory()`, capped on
Linux by the current process cgroup's `free_memory` when that value is
available. Apply an 80% ceiling in two stages: check the packed matrix estimate
from the input count before parsing structures, then check the combined matrix
and selected-coordinate estimate after the first structure reveals the atom
count but before loading the remaining structures. If the platform is
unsupported, the query is zero, or the query fails, skip the hard ceiling and
warn before work starts only when the estimate exceeds 8 GiB.

This is a guard, not a promise that allocation will succeed. Available memory
is a point-in-time estimate that can change immediately. The second estimate
includes only the packed matrix and compact selected-coordinate arrays; it
excludes parser transients, atom-identity validation data, allocator overhead,
output DataFrames, and clustering scratch space. Report the included and
excluded categories, and retain an explicit `bypass_mem_check` override that
skips both warnings and errors. Use available rather than total memory: total
memory includes RAM already committed to other work, while available memory
includes reclaimable capacity and is the closer measure of whether this new
allocation is prudent.

## Accepted design

The completed design interview selected k-medoids only, a fixed internal
automatic minimum of two clusters, an all-zero one-cluster short circuit,
`num_clusters` precedence over `max_clusters`, exact case-sensitive ID sorting,
and `Float64` RMSDs. Two heuristic memory checks and their explicit bypass are
described above. The complete public and operational decision is recorded in
[ADR 0008](../adr/0008-cluster-structures-with-kabsch-and-k-medoids.md).
