# Clustering a pairwise protein-structure RMSD matrix

Research date: 2026-08-28

Alternative methods and crate evidence as of the research date. Accepted choices
are in [ADR 0008](../adr/0008-cluster-structures-with-kabsch-and-k-medoids.md);
arguments, schemas, and measurements are in the
[usage guide](../benchmarks/structure-clustering.md).

## K-medoids and automatic-count evidence

PAM chooses observed medoids to minimize total dissimilarity to the nearest
representative. It does not require a Euclidean feature embedding, and its
observed representatives match the center definition used by GROMACS structure
clustering. FasterPAM improves the SWAP phase by a factor of order `k`; eager
swaps trade exact swap order for comparable quality and more speed.
See [Kaufman and Rousseeuw 1987](https://wis.kuleuven.be/statdatascience/robust/papers/publications-1987/kaufmanrousseeuw-clusteringbymedoids-l1norm-1987.pdf),
[Schubert and Rousseeuw 2021](https://arxiv.org/abs/2008.05171), and
[GROMACS `gmx cluster`](https://manual.gromacs.org/current/onlinehelp/gmx-cluster.html).

`kmedoids` 0.5.5 supports a `LowerTriangle` adapter without an `ndarray` square
copy and is GPL-3.0-or-later. Its authors report parallel FasterPAM as typically
faster above roughly 5,000 observations. That feature also enables `rayon`,
`rand`, and `ndarray`; a future parallel-clustering proposal should include an
end-to-end benchmark above that size. See the
[crate API](https://docs.rs/kmedoids/0.5.5/kmedoids/),
[array adapters](https://docs.rs/kmedoids/0.5.5/kmedoids/arrayadapter/index.html),
[manifest](https://docs.rs/crate/kmedoids/0.5.5/source/Cargo.toml.orig), and
[upstream README](https://github.com/kno10/rust-kmedoids).

DynMSC starts from a maximum medoid count, optimizes average medoid silhouette,
removes medoids in stages, and retains the best solution within the bounds.
This is a criterion-dependent choice, not parameter-free discovery of biological
states. The authors' Python wrapper warns that a high maximum can favor many
singletons and recommends a maximum two to three times the expected count.
See the [DynMSC API](https://docs.rs/kmedoids/0.5.5/kmedoids/fn.dynmsc.html),
[Lenssen and Schubert 2024](https://www.sciencedirect.com/science/article/pii/S0306437923001266),
and [authors' guidance](https://github.com/kno10/python-kmedoids).

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

## Additional Rust crate survey

`linfa-hierarchical` delegates to `kodama`, described above, but adds a
similarity-kernel/`ndarray` interface rather than accepting the existing condensed
RMSD matrix ([Linfa documentation](https://docs.rs/linfa-hierarchical/0.8.1/linfa_hierarchical/)).

No maintained crate directly implementing the Daura/GROMOS conformation
algorithm over a precomputed dissimilarity matrix was found in the crates.io
index under either
[`gromos`](https://crates.io/search?q=gromos) or
[`daura`](https://crates.io/search?q=daura). This is distinct from the
algorithm being unavailable generally: GROMOS++ documents its `cluster`
program as consuming an RMSD matrix, repeatedly choosing the structure with
the most neighbors below a cutoff, and removing that cluster
([GROMOS++ manual](https://www.gromos.net/gromos11_pdf_manuals/vol5.pdf)).
These findings support the k-medoids-only decision in ADR 0008.

## Memory-query alternatives

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
