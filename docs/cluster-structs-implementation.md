# Structure RMSD and clustering implementation plan

This is the implementation ledger for ADR 0008. Work proceeds in vertical,
testable slices and keeps the numerical kernel independent of structure
selection, table I/O, Python, and CLI concerns.

## 1. Dependencies and core types

- [x] Add `kmedoids = { version = "0.5.5", default-features = false }` and a
      direct compatible `sysinfo` dependency; update the lockfile.
- [x] Add Arpeggia-owned atom-selection, structure-observation, packed-matrix,
      clustering-method, and options/result types without exposing dependency
      types.
- [x] Re-export only the public Rust surface from `src/lib.rs`.

## 2. Kabsch RMSD tracer bullet

- [x] Implement one serial `f64` Kabsch kernel over equal-shaped coordinate
      slices with proper-rotation correction and no scaling.
- [x] Reject unequal shapes, non-finite coordinates/results, and fewer than
      three non-collinear points with typed errors; normalize anchor-relative
      displacements with a scaled-subtraction overflow fallback, cache their
      centroids, and use a shared physical magnitude plus a scaled residual norm.
- [x] Test identity, translation, rotation, reflection, noise, large offsets and
      magnitudes, coplanar/collinear inputs, symmetry, and known RMSD values.
- [x] Add the Rust and Python scalar APIs and the `rmsd` CLI command, including
      scalar stdout and ordinary warning/error boundaries.

## 3. Structure selection and correspondence

- [x] Parse the chain/residue union grammar, including full chains, insertion
      codes, bare maximal upper insertion bounds, negative author numbering,
      overlaps, and malformed clauses.
- [x] Implement `ca`, `backbone`, `heavy`, and `all` presets with the approved
      protein/cap/exclusion rules and C-alpha default.
- [x] Apply existing model/conformer selection and diagnostics, construct one
      canonical atom-identity order, reject mismatches against it, and discard
      non-reference identity tables after validation.
- [x] Leave a focused `TODO:` at the correspondence seam for future sequence
      and/or structural alignment; add no speculative alignment abstraction.
- [x] Test multi-chain range selection, terminal caps, hydrogen presence,
      alternate conformers, multiple models, missing atoms, and identity/order
      mismatches.

## 4. Ensemble inventory and tabular input

- [x] Discover case-insensitive PDB/mmCIF extensions non-recursively and derive
      case-sensitive IDs from filename stems.
- [x] Read CSV, Parquet, JSON, and NDJSON manifests with configurable ID/path
      columns, eager CSV/Parquet projection, lazy NDJSON projection pushdown,
      ordinary-JSON projection after parsing, and manifest-relative paths.
- [x] Validate supported extensions, path existence, minimum counts, empty or
      duplicate IDs, and duplicate canonical paths before calculations.
- [x] Sort observations by exact ID and test all input formats, duplicate and
      case behavior, relative paths, and unsupported files.

## 5. Packed pairwise matrix and memory guard

- [x] Implement checked packed-triangle sizing/indexing and strict long-table
      conversion/validation without square copies or temporary ID vectors;
      scan fragmented Polars columns sequentially, pre-size string builders,
      and move packed values when the matrix is consumed.
- [x] Query available/cgroup memory narrowly and apply the matrix-only preflight
      before structure parsing.
- [x] Prepare the first structure, estimate its complete coordinate store, and
      apply the combined second check before preparing the rest.
- [x] Implement `bypass_mem_check`, overflow-safe estimates, the unavailable-RAM
      8-GiB warning, and diagnostics that name included/excluded memory.
- [x] Test boundary arithmetic and guard behavior independently of host RAM by
      injecting memory values into the private decision helper.

## 6. Bounded preparation and pair parallelism

- [x] Prepare the first structure serially, then prepare the rest with
      `min(num_threads, 8)` Rayon workers; retain only normalized-centered
      selected coordinates, scales, and canonical diagnostics, and
      cooperatively stop collection after a parsing failure.
- [x] Fill disjoint contiguous packed-matrix chunks with at most one worker per
      pair, keeping each Kabsch calculation serial and lock-free; cache each
      structure's numerical coordinate scale outside the pair loop.
- [x] Expose `pairwise_rmsd` through Rust and Python and return the three-column
      long Polars DataFrame.
- [x] Verify identical matrix values and ordering across thread counts and
      deterministic warning order across parallel preparation.

## 7. K-medoids clustering

- [x] Adapt the packed matrix to `kmedoids` without `ndarray` or a square copy.
- [x] Apply only the uniform divisor needed to keep PAM aggregate losses finite,
      reject distance ranges that collapse relative to their maximum, and
      retain original Angstrom values for output; avoid arithmetic in ordinary
      unscaled adapter lookups.
- [x] Implement PAM BUILD plus FasterPAM for fixed counts and DynMSC for bounded
      automatic counts, including argument precedence, numerically identical
      ensembles, and verified iteration failures.
- [x] Normalize labels deterministically, resolve fixed-count equal-loss medoid
      ties by canonical input order with reoptimization, preserve DynMSC's
      silhouette objective, and short-circuit all-zero matrices to one cluster.
- [x] Produce only `id`, `cluster_id`, `medoid_id`, and `rmsd_to_medoid` with
      the approved Polars dtypes.
- [x] Test fixed `k=1`, `k=n`, automatic bounds, identical structures, tied
      and duplicate structures, row permutation, invalid distances, extreme
      finite distance scales, expected
      automatic partitions and medoids, and iteration limits. Automatic bounds
      exclude `n` to avoid a trivially optimal all-singleton silhouette.

## 8. CLI output and ID-only reuse

- [x] Add `cluster-structs` with directory input, selection/model,
      cluster bound, iteration, output, threading, and memory-bypass arguments.
- [x] Validate count-independent cluster bounds and residue syntax before CLI
      directory discovery or Python structure loading.
- [x] Write optional `id_1`, `id_2`, `rmsd` output before clustering in CSV,
      Parquet, JSON, or NDJSON, then release any transient output DataFrame.
- [x] Check only the exact requested pairwise path and reuse a complete table
      with the exact current ID set; emit a debug log on reuse.
- [x] Reject wrong-size cached Parquet tables from row-count metadata before
      decoding and bound cached CSV, JSON, and NDJSON reads to one row beyond
      the expected count.
- [x] Fail without overwriting when an existing cache is malformed, incomplete,
      or ID-mismatched, and instruct the caller to remove it.
- [x] Test successful recovery/reuse, stale-selection reuse as the documented
      ID-only behavior, forced recalculation by removal, output ordering, and a
      clustering failure after successful pairwise persistence.

## 9. Python and documentation contract

- [x] Add Python `rmsd`, `pairwise_rmsd`, and `cluster_structs`; require exactly
      one of structure input or a Polars pairwise DataFrame for clustering.
- [x] Reject count-independent clustering options and invalid RMSD atom presets
      before structure parsing or pairwise-table packing.
- [x] Update Python exports, contracts, type stubs, docstrings, and exception and
      warning tests.
- [x] Document CLI examples, selection grammar, formats, output schemas,
      quadratic scaling, both memory checks and their exclusions, bounded
      parser/RMSD threading, single-threaded clustering, and ID-only cache risk.
- [x] Update the changelog only when the public implementation is complete.

## 10. Performance and completion gate

- [x] Benchmark preparation with 1 versus `min(num_threads, 8)` workers and
      pairwise construction with 1, 2, 4, 8, and automatic workers across
      representative structure and atom counts.
- [x] Record end-to-end time, pairwise time, clustering time, and peak RSS;
      verify that output serialization is measured separately.
- [x] Keep k-medoids serial unless measurements identify it as a serious
      end-to-end bottleneck; do not enable its parallel feature speculatively.
- [x] Run `cargo fmt`, Clippy, Rust tests, doctests, CLI integration tests,
      Python tests, `ty`, and `prek` on every changed file.
- [x] Audit the implementation against ADR 0008 and this ledger before marking
      the feature complete.
