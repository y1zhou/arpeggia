# Structure RMSD and clustering

Arpeggia superposes exactly corresponding protein atoms with the Kabsch
algorithm and clusters the resulting pairwise RMSD matrix with k-medoids. This
release does not align sequences or infer missing-atom correspondence: selected
chain IDs, residue identities, and atom names must match exactly.

## RMSD

The default selection is every coordinate-observed protein residue and its
C-alpha atom. `atoms` accepts `ca`, `backbone` (`N`, `CA`, `C`, `O`, and
`OXT`), `heavy`, or `all`. Heavy selection excludes hydrogen and deuterium,
including digit-leading atom names when element metadata is absent. Heavy and
all-atom selections retain ACE/NH2 caps but exclude solvent, ions, and ligands.

Residue selection is a comma-separated union of chain and author-residue
clauses. A bare chain selects all its residues. For example,
`A:1-100,A:110-120,B,C:1,C:3,C:5,C:7,C:9-20` excludes A:101-109, includes all
of B, and selects the listed parts of C. Negative numbering and insertion codes
are valid, such as `A:-5--1` and `B:10A-20`.
A bare upper bound includes every insertion code at that author residue, so
`A:10A-10` selects insertion 10A through the final insertion at residue 10.

```bash
arpeggia rmsd reference.cif mobile.cif \
  --residues "A:1-120,B" --atoms backbone
```

The CLI prints one RMSD in Ångströms. Python provides the scalar operation and
the complete unordered pair table:

```python
import arpeggia

value = arpeggia.rmsd("reference.cif", "mobile.cif", atoms="ca")
pairs = arpeggia.pairwise_rmsd(
    "structures/", residues="A:1-120,B", atoms="ca", num_threads=8
)
```

`pairwise_rmsd` accepts a non-recursive structure directory or a CSV, Parquet,
JSON, or NDJSON manifest. Manifest columns default to `id` and `path` and can
be changed with `id_col` and `path_col`; relative paths resolve against the
manifest. Directory IDs are case-sensitive filename stems, while PDB/mmCIF
extensions are case-insensitive. Duplicate IDs or canonical paths fail before
RMSD calculation. Structure, manifest, and cache paths must be regular files.
The result has one unordered pair per row:

| column | type |
| --- | --- |
| `id_1` | String |
| `id_2` | String |
| `rmsd` | Float64 |

## Clustering

Use either a fixed cluster count or an automatically selected count bounded by
`max_clusters`. A fixed count uses deterministic PAM BUILD initialization and
FasterPAM. Automatic selection uses DynMSC over 2 through `max_clusters`; an
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

The CLI intentionally accepts only a non-recursive structure directory. Python
accepts exactly one of a directory/manifest `input` or a complete long-form
Polars `pairwise_rmsd` DataFrame, allowing a calculated matrix to be reused
without recomputation:

```python
clusters = arpeggia.cluster_structs(
    pairwise_rmsd=pairs,
    max_clusters=10,
)
```

The cluster table contains:

| column | type | meaning |
| --- | --- | --- |
| `id` | String | input structure ID |
| `cluster_id` | UInt32 | deterministic zero-based cluster label |
| `medoid_id` | String | observed representative structure |
| `rmsd_to_medoid` | Float64 | RMSD to that representative |

CLI tables can be CSV, Parquet, JSON, or NDJSON. `--pairwise-rmsd` writes the
pair table before clustering, preserving it if clustering fails. A later run
reuses the exact requested pairwise path only when its schema, complete pair
coverage, and ID set validate. Cache reuse checks IDs only—not file contents,
selection, model, conformers, or Arpeggia version. Remove the pairwise file to
force recalculation. Wrong-size caches are rejected with bounded reads before
their complete tables are materialized.

## Memory and threads

Pairwise RMSD is quadratic in structure count. Before parsing coordinates,
Arpeggia estimates the packed matrix as `4n(n-1)` bytes. After preparing the
first structure, it also estimates selected coordinates as `24na` bytes for
`n` structures and `a` selected atoms. Either estimate fails above 80% of
effective available RAM; `bypass_mem_check=True` or `--bypass-mem-check`
disables this heuristic.

These estimates cover only the packed RMSD matrix and selected coordinate
arrays. They exclude full-structure parser transients, atom-identity keys,
allocator overhead, output DataFrames, serialization, and clustering scratch
space, so they are not a maximum-RAM guarantee.

The first structure is prepared serially. Remaining structures use at most
`min(num_threads, 8)` parser workers to avoid saturating storage; pairwise RMSD
uses up to the smaller of the requested worker count and number of pairs. Each
Kabsch solve and k-medoids clustering remains single-threaded. `num_threads=0`
selects available processors.

The implementation decision and local performance measurements are recorded in
[ADR 0008](adr/0008-cluster-structures-with-kabsch-and-k-medoids.md) and the
[local benchmark report](benchmarks/cluster-structs-local.md).
