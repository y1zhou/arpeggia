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
| heavy | 1 | 1.10 | 148.8 |
| heavy | 8 | 0.18 | 173.4 |

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
eight. Peak RSS substantially exceeds the documented packed-matrix plus
coordinate estimate for heavy atoms because the guard intentionally excludes
parser transients, identity keys, allocator overhead, DataFrames, and
clustering scratch space.
