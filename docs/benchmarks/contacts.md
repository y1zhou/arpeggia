# Contact analysis benchmarks

Historical measurements from the v0.9.2 cleanup. The NDJSON reader and
ring-ring enumeration were measured separately; neither comparison establishes
a speedup for every stage or input. See also the
[package-size measurements](https://github.com/y1zhou/arpeggia/blob/master/docs/benchmarks/package-size.md).

## 5B8C contacts: lazy versus eager NDJSON

### Result

The eager replacement showed no material slowdown on this input. Contact
generation differed by 0.8% with one thread and 1.2% with eight threads.
One-thread NDJSON reads were effectively unchanged when checking row counts
and 9–13% faster without that check. With eight threads, eager reads were
51–65% faster. This supports retaining the removal of the lazy engine.

### Timings

Milliseconds below are the median of seven batch medians. Negative change
means the eager implementation was faster.

| Threads | Operation | Lazy (ms) | Eager (ms) | Change |
| ---: | --- | ---: | ---: | ---: |
| 1 | Generate contacts and serialize NDJSON | 24.038 | 24.220 | +0.8% |
| 1 | Read all 20 columns | 2.286 | 2.073 | -9.3% |
| 1 | Read all 20 columns + row-count check | 2.267 | 2.278 | +0.5% |
| 1 | Read 3 columns | 1.942 | 1.694 | -12.8% |
| 1 | Read 3 columns + row-count check | 1.926 | 1.915 | -0.6% |
| 8 | Generate contacts and serialize NDJSON | 19.726 | 19.972 | +1.2% |
| 8 | Read all 20 columns | 1.500 | 0.593 | -60.4% |
| 8 | Read all 20 columns + row-count check | 1.456 | 0.711 | -51.2% |
| 8 | Read 3 columns | 1.229 | 0.432 | -64.9% |
| 8 | Read 3 columns + row-count check | 1.238 | 0.592 | -52.2% |

### Input and environment

- Input: unhydrogenated `5B8C.pdb`, 469,557 bytes.
- Input SHA-256: `574a7ad1a4846aaae3f8e2924d67ce7a34fa76b0e0acf8fee9471ec1e4a56766`.
- Generated contacts: 2,574 rows × 20 columns; 969,835 bytes of NDJSON.
- Contact options: groups `/`, VdW compensation 0.1 Å, distance cutoff 6.5 Å, `AllCharged` protonation, pH 7.4.
- Host: AMD Ryzen 9 9950X3D, 16 cores / 32 logical CPUs, Linux x86-64.
- Rust: `rustc 1.96.0 (ac68faa20 2026-05-25)`; Polars 0.55.2; optimized release builds.
- Both `POLARS_MAX_THREADS` and `RAYON_NUM_THREADS` were set to the reported thread count. Contact analysis used the same explicit count.

### Method

Two binaries were built from the same v0.9.2 cleanup source snapshot. The eager
variant used the cleanup implementation. The lazy variant enabled
`polars/lazy` and used the original `read_dataframe` function from commit
`7e06c82549c1eb94a1e75b2e1d5c2d38f426bc99`. Scientific code and all other
dependency versions were identical. This isolates the feature/reader change
from the other cleanup edits.

Each operation ran in seven batches per variant, with variant order shuffled
using seed 5808. Each subprocess performed an untimed validation call and
three warmups, then 15 timed iterations for contacts or 50 for reads.
That gives 105 contact samples or 350 reader samples per variant/thread
configuration. Timing used Rust `Instant` around the operation; subprocess
startup and result destruction were outside the timer.

Contact timing includes PDB loading, metadata loading, contact calculation,
and NDJSON serialization to an I/O sink. It excludes CLI logging and disk
output. Reader timing includes the production file checks, schema inference,
projection, and materialization. The projected columns were `from_chain`,
`to_chain`, and `distance`. Bounded reads supplied the actual expected row
count, 2,574, exercising the eager preflight scan and the lazy row limit.

Every batch checked the result shape, schema, and checksum of serialized
NDJSON. They matched between variants for every operation and thread count.
The contacts table serves as a realistic NDJSON input; contact generation
itself does not call the NDJSON reader in production.

### Limits

These are steady-state measurements with warm filesystem caches on a small
local table. They do not establish performance for cold disks, remote
filesystems, or much larger RMSD caches. The eager row-count scan remains an
additional sequential pass. The near-1% generation differences are small
absolute differences (0.18–0.25 ms), not evidence that the scientific
algorithm became slower.

The benchmark harness and raw samples were not committed; the summary timings
and methodology above are the retained evidence.

## Ring-ring enumeration benchmark

This compares buffered candidate indices with borrowed pair enumeration inside
the parallel ring loop: both variants include shared plane
preparation and all other follow-up changes. Timed work starts after preparation
and includes ring-ring classification and output disposal, not parsing, plane
fitting, atom contacts, or DataFrame construction. Runs used Rust 1.96.0 release
builds on Linux x86_64, AMD Ryzen 9 9950X3D, with Rayon limited to one or eight
threads.

The stress input contains 128 translated copies of the 44 prepared 5B8C rings,
with distinct residue identities and corresponding peptide-neighbor exclusions.
Copies are separated by 200 Å. This is a synthetic ring-only scaling case, not a
larger experimentally observed protein or a benchmark of total contact analysis.
No rings or atoms are reconstructed inside the timed pair loop.

| Input | Rings | Threads | Buffered pairs (ms) | Streamed pairs (ms) | Change | Peak RSS before/after (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 5B8C | 44 | 1 | 0.03712 | 0.02115 | −43.0% | 4.96 / 4.95 |
| 5B8C | 44 | 8 | 0.03744 | 0.00755 | −79.8% | 5.03 / 4.95 |
| 6BFT | 138 | 1 | 0.21659 | 0.16127 | −25.5% | 8.43 / 8.48 |
| 6BFT | 138 | 8 | 0.20211 | 0.03424 | −83.1% | 8.43 / 8.22 |
| 5B8C × 128, synthetic | 5,632 | 1 | 514.19612 | 421.44474 | −18.0% | 249.02 / 7.23 |
| 5B8C × 128, synthetic | 5,632 | 8 | 467.48880 | 55.31827 | −88.2% | 248.73 / 6.99 |

Times are medians of seven process medians, alternating variant order. Each
process performs an untimed warmup, then 1,000 iterations for real structures or
three for the stress case. Peak RSS is the median process high-water mark from
GNU `time`; it includes setup and is not an allocation-only measurement.
The variants have identical sorted-output fingerprints and row counts: 8 for
5B8C, 43 for 6BFT, and 1,024 for the synthetic case. Small-input RSS differences
are negligible; the synthetic case demonstrates the removed quadratic storage.
Pair enumeration itself remains quadratic. These results support retaining the
ring-ring change without extending it to the other candidate paths.

The isolated harness and raw measurements were not committed.

Across the complete follow-up cleanup, contact CSVs matched baseline
`4481556` byte for byte for 5B8C (2,574 rows), hydrogenated 5B8C (2,645),
and 6BFT with all-to-all (7,181), `H,L/C,G` (128), and overlapping
`H,L/H,C` (1,629) selections. This checks combined output preservation;
the timings above isolate only pair enumeration.
