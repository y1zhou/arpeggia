# Structure RMSD and clustering implementation plan

This is the implementation ledger for ADR 0008. Work proceeds in vertical,
testable slices and keeps the numerical kernel independent of structure
selection, table I/O, Python, and CLI concerns.

## 1. Dependencies and core types

- [ ] Add `kmedoids = { version = "0.5.5", default-features = false }` and a
      direct compatible `sysinfo` dependency; update the lockfile.
- [ ] Add Arpeggia-owned atom-selection, structure-observation, packed-matrix,
      clustering-method, and options/result types without exposing dependency
      types.
- [ ] Re-export only the public Rust surface from `src/lib.rs`.

## 2. Kabsch RMSD tracer bullet

- [ ] Implement one serial `f64` Kabsch kernel over equal-shaped coordinate
      slices with proper-rotation correction and no scaling.
- [ ] Reject unequal shapes, non-finite coordinates, and fewer than three
      non-collinear points with typed errors.
- [ ] Test identity, translation, rotation, reflection, noise, large offsets,
      coplanar/collinear inputs, symmetry, and known RMSD values.
- [ ] Add the Rust and Python scalar APIs and the `rmsd` CLI command, including
      scalar stdout and ordinary warning/error boundaries.

## 3. Structure selection and correspondence

- [ ] Parse the chain/residue union grammar, including full chains, insertion
      codes, negative author numbering, overlaps, and malformed clauses.
- [ ] Implement `ca`, `backbone`, `heavy`, and `all` presets with the approved
      protein/cap/exclusion rules and C-alpha default.
- [ ] Apply existing model/conformer selection and diagnostics, construct one
      canonical atom-identity order, and reject mismatches against it.
- [ ] Leave a focused `TODO:` at the correspondence seam for future sequence
      and/or structural alignment; add no speculative alignment abstraction.
- [ ] Test multi-chain range selection, terminal caps, hydrogen presence,
      alternate conformers, multiple models, missing atoms, and identity/order
      mismatches.

## 4. Ensemble inventory and tabular input

- [ ] Discover case-insensitive PDB/mmCIF extensions non-recursively and derive
      case-sensitive IDs from filename stems.
- [ ] Read CSV, Parquet, JSON, and NDJSON manifests with configurable ID/path
      columns and manifest-relative paths.
- [ ] Validate supported extensions, path existence, minimum counts, empty or
      duplicate IDs, and duplicate canonical paths before calculations.
- [ ] Sort observations by exact ID and test all input formats, duplicate and
      case behavior, relative paths, and unsupported files.

## 5. Packed pairwise matrix and memory guard

- [ ] Implement checked packed-triangle sizing/indexing and strict long-table
      conversion/validation without square or duplicate matrix copies.
- [ ] Query available/cgroup memory narrowly and apply the matrix-only preflight
      before structure parsing.
- [ ] Prepare the first structure, estimate its complete coordinate store, and
      apply the combined second check before preparing the rest.
- [ ] Implement `bypass_mem_check`, overflow-safe estimates, the unavailable-RAM
      8-GiB warning, and diagnostics that name included/excluded memory.
- [ ] Test boundary arithmetic and guard behavior independently of host RAM by
      injecting memory values into the private decision helper.

## 6. Bounded preparation and pair parallelism

- [ ] Prepare the first structure serially, then prepare the rest with
      `min(num_threads, 8)` Rayon workers; retain only centered selected
      coordinates and canonical diagnostics.
- [ ] Fill disjoint contiguous packed-matrix chunks with all requested workers,
      keeping each Kabsch calculation serial and lock-free.
- [ ] Expose `pairwise_rmsd` through Rust and Python and return the three-column
      long Polars DataFrame.
- [ ] Verify identical matrix values and ordering across thread counts and
      deterministic error/warning order across parallel preparation.

## 7. K-medoids clustering

- [ ] Adapt the packed matrix to `kmedoids` without `ndarray` or a square copy.
- [ ] Implement PAM BUILD plus FasterPAM for fixed counts and DynMSC for bounded
      automatic counts, including argument precedence and iteration failures.
- [ ] Normalize labels and medoid ties deterministically and short-circuit
      all-zero matrices to one cluster.
- [ ] Produce only `id`, `cluster_id`, `medoid_id`, and `rmsd_to_medoid` with
      the approved Polars dtypes.
- [ ] Test fixed `k=1`, `k=n`, automatic bounds, identical structures, tied
      distances, row permutation, invalid distances, and iteration limits.

## 8. CLI output and ID-only reuse

- [ ] Add `cluster-structs` with directory-or-manifest input, selection/model,
      cluster bound, iteration, output, threading, and memory-bypass arguments.
- [ ] Write optional `id_1`, `id_2`, `rmsd` output before clustering in CSV,
      Parquet, JSON, or NDJSON, then release any transient output DataFrame.
- [ ] Check only the exact requested pairwise path and reuse a complete table
      with the exact current ID set; emit a debug log on reuse.
- [ ] Fail without overwriting when an existing cache is malformed, incomplete,
      or ID-mismatched, and instruct the caller to remove it.
- [ ] Test successful recovery/reuse, stale-selection reuse as the documented
      ID-only behavior, forced recalculation by removal, output ordering, and a
      clustering failure after successful pairwise persistence.

## 9. Python and documentation contract

- [ ] Add Python `rmsd`, `pairwise_rmsd`, and `cluster_structs`; require exactly
      one of structure input or a Polars pairwise DataFrame for clustering.
- [ ] Update Python exports, contracts, type stubs, docstrings, and exception and
      warning tests.
- [ ] Document CLI examples, selection grammar, formats, output schemas,
      quadratic scaling, both memory checks and their exclusions, bounded
      parser/RMSD threading, single-threaded clustering, and ID-only cache risk.
- [ ] Update the changelog only when the public implementation is complete.

## 10. Performance and completion gate

- [ ] Benchmark preparation with 1 versus `min(num_threads, 8)` workers and
      pairwise construction with 1, 2, 4, 8, and automatic workers across
      representative structure and atom counts.
- [ ] Record end-to-end time, pairwise time, clustering time, and peak RSS;
      verify that output serialization is measured separately.
- [ ] Keep k-medoids serial unless measurements identify it as a serious
      end-to-end bottleneck; do not enable its parallel feature speculatively.
- [ ] Run `cargo fmt`, Clippy, Rust tests, doctests, CLI integration tests,
      Python tests, `ty`, and `prek` on every changed file.
- [ ] Audit the implementation against ADR 0008 and this ledger before marking
      the feature complete.
