# Local structure-clustering benchmark

This local, doc-only benchmark used the 250 mmCIF files in
`/tmp/arpeggia-benchmark-structs/` (52 MiB) on 2026-08-28. The host exposed 32
processors and 123 GiB RAM. Commands used a locked release build, fixed
`k=5`, and three warm-cache repetitions in shuffled order. `/usr/bin/time`
reported end-to-end wall time and process peak RSS; benchmark outputs remained
under `/tmp` and are not packaged or committed.

## Default C-alpha selection

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

## Atom-count and stage proxies

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

## Interpretation

The bounded eight-worker preparation path and parallel pair filling materially
reduce runtime. Eight workers are the practical default for this dataset:
automatic 32-worker pair calculation adds no speed but has similar memory to
eight. Retaining atom identities only for the reference structure and releasing
parser workers before pairwise calculation reduced the heavy-atom peak RSS
from 148.8 to 35.5 MiB with one worker and from 173.4 to 57.0 MiB with eight
workers, without changing either output. Peak RSS still
exceeds the documented packed-matrix plus coordinate estimate because the
guard intentionally excludes parser transients, the reference identity keys,
allocator overhead, DataFrames, and clustering scratch space.

## Independent superposition and RMSD selections

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
