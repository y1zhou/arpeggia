# Structure RMSD and clustering

Arpeggia superposes exactly corresponding protein atoms with the Kabsch
algorithm and clusters the resulting pairwise RMSD matrix with k-medoids. It
does not align sequences or infer missing-atom correspondence: selected
chain IDs, residue identities, and atom names must match exactly.

## RMSD

The Superposition Selection determines the Kabsch transform. The RMSD Selection
is evaluated after applying that fixed transform and does not influence the fit.
Both independently default to every coordinate-observed amino acid recognized by
Arpeggia, so callers pass the same subset to both arguments when they want the
traditional fit-and-score-the-same-atoms calculation.

One shared `atoms` argument accepts `ca`, `backbone` (`N`, `CA`, `C`, `O`, and
`OXT`), `heavy`, or `all`. Heavy selection excludes hydrogen and deuterium,
including digit-leading atom names when element metadata is absent. Heavy and
all-atom selections retain ACE/NH2 caps but exclude solvent, ions, and ligands.
Arbitrary polymers, ligands, and modified residues are not yet retained even in
`all` mode.

Residue selection is a comma-separated union of chain and author-residue
clauses. A bare chain selects all its residues. For example,
`A:1-100,A:110-120,B,C:1,C:3,C:5,C:7,C:9-20` excludes A:101-109, includes all
of B, and selects the listed parts of C. Negative numbering and insertion codes
are valid, such as `A:-5--1` and `B:10A-20`.
A bare upper bound includes every insertion code at that author residue, so
`A:10A-10` selects insertion 10A through the final insertion at residue 10.

```bash
arpeggia rmsd reference.cif mobile.cif \
  --superpose-residues "A" \
  --rmsd-residues "B,C" \
  --atoms backbone
```

This example establishes the coordinate frame from chain A and reports the
motion of chains B and C relative to it. Superposition requires at least three
non-collinear atom pairs; RMSD evaluation requires at least one atom pair. The
CLI prints one RMSD in Ångströms. Python provides the scalar operation and the
complete unordered pair table:

```python
import arpeggia

value = arpeggia.rmsd(
    "reference.cif",
    "mobile.cif",
    superpose_residues="A",
    rmsd_residues="B,C",
    atoms="ca",
)
pairs = arpeggia.pairwise_rmsd(
    "structures/",
    superpose_residues="A",
    rmsd_residues="B,C",
    atoms="ca",
    num_threads=8,
)
```

`pairwise_rmsd` accepts a non-recursive structure directory or a CSV, Parquet,
or NDJSON manifest. Manifest columns default to `id` and `path` and can
be changed with `id_col` and `path_col`; relative paths resolve against the
manifest. Directory IDs are case-sensitive filename stems, while PDB/mmCIF
extensions are case-insensitive. Duplicate IDs or canonical paths fail before
RMSD calculation. Structure, manifest, and cache paths must be regular files.
Pairwise calculation requires at least two structures. The result has one
unordered pair per row:

| column | type |
| --- | --- |
| `id_1` | String |
| `id_2` | String |
| `rmsd` | Float64 |

## Clustering

Clustering requires at least three structures. Use either a fixed cluster
count in `1..=n` or an automatically selected count bounded by `max_clusters`.
A fixed count uses deterministic PAM BUILD initialization and FasterPAM. Automatic selection uses DynMSC over 2 through `max_clusters`; an
ensemble whose pairwise RMSDs are all at most `1e-12` Angstrom becomes one
deterministic cluster. `max_clusters` must be smaller than the number of
structures because an all-singleton partition has a trivially maximal medoid
silhouette. If both bounds are supplied, the fixed count wins with a warning.
Fixed-count equal-loss medoid ties use canonical input order and are
reoptimized after a tie move. Automatic clustering preserves DynMSC's
medoid-silhouette objective and deterministic canonical input order.

```bash
arpeggia cluster-structs \
  --input structures/ \
  --output results/ \
  --num-clusters 5 \
  --pairwise-rmsd \
  --num-threads 8
```

The CLI accepts only a non-recursive structure directory. Python
accepts exactly one of a directory/manifest `input` or a complete long-form
Polars `pairwise_rmsd` DataFrame, allowing a calculated matrix to be reused
without recomputation:

```python
clusters = arpeggia.cluster_structs(
    pairwise_rmsd=pairs,
    max_clusters=10,
)
```

`max_iterations` defaults to 100. Failure to converge within the supported
iteration budget raises a calculation error.

The cluster table contains:

| column | type | meaning |
| --- | --- | --- |
| `id` | String | input structure ID |
| `cluster_id` | UInt32 | deterministic zero-based cluster label |
| `medoid_id` | String | observed representative structure |
| `rmsd_to_medoid` | Float64 | RMSD to that representative |

CLI tables can be CSV, Parquet, or NDJSON. `--pairwise-rmsd` writes the
pair table before clustering, preserving it if clustering fails. A later run
reuses the exact requested pairwise path only when its schema, complete pair
coverage, and ID set validate. Cache reuse checks IDs only—not file contents,
either residue selection, model, conformers, or Arpeggia version. Remove the
pairwise file to force recalculation. Malformed or ID-mismatched caches fail
without being overwritten; wrong-size caches are rejected before their
complete tables are materialized.

## Memory and threads

Pairwise RMSD is quadratic in structure count. Before parsing coordinates,
Arpeggia estimates the packed matrix as `4n(n-1)` bytes. After preparing the
first structure, it also estimates selected coordinates as `24nu` bytes for
`n` structures and `u` atoms in the union of both selections. Overlapping atoms
are stored once. Either estimate fails above 80% of
effective available RAM; `bypass_mem_check=True` or `--bypass-mem-check`
disables this heuristic. If available RAM cannot be queried, estimates above
8 GiB produce a warning instead of a hard limit.

These estimates cover only the packed RMSD matrix and selected coordinate
arrays. They exclude full-structure parser transients, atom-identity keys,
allocator overhead, output DataFrames, serialization, and clustering scratch
space, so they are not a maximum-RAM guarantee.

The first structure is prepared serially. Remaining structures use at most
`min(num_threads, 8)` parser workers to avoid saturating storage; pairwise RMSD
uses up to the smallest of the requested worker count, available processors,
and number of pairs. Each Kabsch solve and k-medoids clustering remains single-threaded. `num_threads=0`
selects available processors.

Algorithm choices and their rationale are recorded in
[ADR 0008](../adr/0008-cluster-structures-with-kabsch-and-k-medoids.md).

## Local structure-clustering benchmark

Measured on 2026-08-28 using 250 mmCIF files (52 MiB), 32 processors, and
123 GiB RAM. Locked release builds used fixed `k=5` and three warm-cache
repetitions in shuffled order. GNU `time` measured end-to-end wall time and
process peak RSS. The input corpus and generated outputs are not committed.

### Default C-alpha selection

Each run discovered and prepared all 250 structures, calculated 31,125 RMSD
pairs, wrote the 2.9 MiB long pair table, clustered it, and wrote the cluster
table.

| RMSD workers | median time (s) | median peak RSS (MiB) |
| ---: | ---: | ---: |
| 1 | 0.75 | 34.8 |
| 2 | 0.38 | 41.8 |
| 4 | 0.19 | 45.7 |
| 8 | 0.11 | 61.5 |
| automatic (32) | 0.11 | 62.7 |

Eight workers were 6.8× faster than one. Automatic selection did not improve
on eight workers at this problem size. The first cold single-worker run took
2.05 s and was excluded from the warm-cache medians.

The CSV pair tables and cluster tables were byte-identical across 1, 2, 4, 8,
and automatic workers. Pairwise RMSD ranged from 0.172710 to 1.276763 Å, with a
0.481133 Å median. Fixed `k=5` produced cluster sizes 27, 43, 66, 55, and 59.
Automatic DynMSC with `max_clusters=10` selected two clusters and completed the
cache-backed workflow in 0.05 s.

### Atom-count and stage proxies

Heavy-atom selection increased coordinate storage and atom-loop work:

| selection | workers | median time (s) | median peak RSS (MiB) |
| --- | ---: | ---: | ---: |
| heavy | 1 | 1.07 | 35.5 |
| heavy | 8 | 0.15 | 57.0 |

Heavy-atom cluster outputs were byte-identical between worker counts. A
three-C-alpha selection, used as a preparation-dominated proxy while retaining
the public pairwise workflow, took 0.64 s and 31.9 MiB with one worker versus
0.09 s and 43.6 MiB with eight. It is not a pure parser microbenchmark because
the fixed number of Kabsch solves and clustering still run.

Reusing the C-alpha CSV pair table took a 0.04 s median and 31.9 MiB peak RSS
across five runs; this includes directory inventory, pair-table parsing,
k-medoids, and cluster serialization. Serializing the already materialized
31,125-row pair table to CSV with Polars took a 1.19 ms median over ten writes,
so serialization was negligible in the first-run timing.

### Interpretation

The bounded eight-worker preparation path and parallel pair filling materially
reduce runtime. Eight workers are the practical default for this dataset:
automatic 32-worker pair calculation adds no speed but has similar memory to
eight. Retaining atom identities only for the reference structure and releasing
parser workers before pairwise calculation reduced the heavy-atom peak RSS
from 148.8 to 35.5 MiB with one worker and from 173.4 to 57.0 MiB with eight
workers, without changing either output. Peak RSS exceeds the estimate under
[Memory and threads](#memory-and-threads).

### Independent superposition and RMSD selections

The independent-selection implementation was checked on the same 250
structures on 2026-09-01. A preserved pre-change release binary and the new
release binary were stripped identically before measurement so debug-section
layout did not distort process RSS. Each median below covers five warm-cache
runs with fixed `k=5`; pair tables were not serialized during timing.

For equal all-residue selections, the pre-change and new pairwise RMSD CSVs and
cluster CSVs were byte-identical. The generalized fixed-transform path also
matches the existing prepared Kabsch path bit-for-bit when its parsed selectors
are equal, while production dispatch retains the existing fast path.

| implementation | workers | median time (s) | median peak RSS (MiB) |
| --- | ---: | ---: | ---: |
| pre-change, equal | 1 | 1.38 | 31.5 |
| new, equal | 1 | 1.37 | 30.8 |
| pre-change, equal | 8 | 0.19 | 46.2 |
| new, equal | 8 | 0.19 | 45.2 |

The new equal-selection path changed median runtime by -0.7% with one worker
and 0.0% with eight, within the 5% and 10% gates. Median peak RSS changed by
-2.4% and -2.1%, within the 10% gate.

The selected C-alpha counts were 193 for chain A and 119 for chain H. The
overlap case fitted `A:1-193,H:1-60` and evaluated
`A:50-193,H:1-119`; the disjoint case fitted A and evaluated H.

| selection | fit `f` | score `r` | union `u` | workers | median time (s) | median peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| equal | 312 | 312 | 312 | 1 | 1.37 | 30.8 |
| equal | 312 | 312 | 312 | 8 | 0.19 | 45.2 |
| overlapping | 253 | 263 | 312 | 1 | 1.51 | 30.9 |
| overlapping | 253 | 263 | 312 | 8 | 0.21 | 45.4 |
| disjoint | 193 | 119 | 312 | 1 | 1.22 | 31.1 |
| disjoint | 193 | 119 | 312 | 8 | 0.17 | 45.5 |

All three cases retain exactly `24nu = 1,872,000` coordinate bytes. Compared
with storing fit and score arrays independently, union storage saves
`24n(f+r-u)`: 1,872,000 bytes for equal selections, 1,224,000 bytes for the
overlapping selections, and zero for disjoint selections. These payload figures
cover only coordinates; observed RSS also includes structure parsing, atom
identities, the packed pair matrix, allocation overhead, and clustering.
