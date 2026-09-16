# Python and CLI examples

Replace file paths and chain IDs with those in your structures. Python examples
use `import arpeggia`. CLI table commands take an output directory with `-o` and
support `-t csv`, `-t parquet` or `-t ndjson`.

## Contacts and table output

Compare chains A/B against C/D, omit zero-occupancy atoms, and choose histidine
typing explicitly. `ph` affects only `protonation="heuristic"`; alternatives are
`"all-charged"` (default) and `"explicit-only"`. Distances are in Å.

```python
import arpeggia

contacts = arpeggia.contacts(
    "structure.pdb", groups="A,B/C,D", ignore_zero_occupancy=True,
    vdw_comp=0.1, dist_cutoff=6.5, protonation="heuristic", ph=7.4,
)
contacts.write_csv("contacts.csv")
contacts.write_parquet("contacts.parquet")
contacts.write_ndjson("contacts.ndjson")
```

```bash
arpeggia contacts -i structure.pdb -o results/ -g "A,B/C,D" \
  --ignore-zero-occupancy --vdw-comp 0.1 --dist-cutoff 6.5 \
  --protonation heuristic --ph 7.4 -t parquet
```

Use `groups="A/"` for A against every other chain, or `"/"` for all inter- and
intra-chain contacts. Results are Polars DataFrames; optional pandas conversion
uses `contacts.to_pandas()` and requires pandas and PyArrow. See
[contact-table recipes](https://github.com/y1zhou/arpeggia/blob/master/docs/scientific-conventions.md#contact-table-examples)
for hydrogen-bond counts and interface residues.

## Surface measurements

SASA levels are `atom`, `residue` and `chain`. `chains` is a comma-separated
string; an empty string includes every chain. `model_num=0` / `--model 0` selects
the first model; other values select model serials.

```python
for level in ("atom", "residue", "chain"):
    surface = arpeggia.sasa(
        "structure.pdb", level=level, chains="H,L", model_num=0,
        probe_radius=1.4, n_points=100,
    )
    surface.write_csv(f"sasa_{level}.csv")

rsa = arpeggia.relative_sasa("structure.pdb", chains="H,L", probe_radius=1.4)
sap = arpeggia.sap_score(
    "structure.pdb", level="residue", chains="H,L",
    probe_radius=1.1, sap_radius=5.0,
)
print(sap.sort("sap_score", descending=True).head(10))
```

```bash
for level in atom residue chain; do
  arpeggia sasa -i structure.pdb -o results/ --filename "sasa_$level" \
    --level "$level" --chains H,L --model 0 --probe-radius 1.4 --num-points 100
done
arpeggia relative-sasa -i structure.pdb -o results/ --chains H,L
arpeggia sap -i structure.pdb -o results/ --level residue --chains H,L \
  --probe-radius 1.1 --sap-radius 5.0
```

Probe and SAP neighborhood radii are in Å. Relative SASA normalizes residue
SASA by reference maxima; SAP describes hydrophobic exposure in a spatial
neighborhood. Their assumptions are in
[Scientific Conventions](https://github.com/y1zhou/arpeggia/blob/master/docs/scientific-conventions.md).

## Buried area and shape complementarity

These calculations require two disjoint, non-empty chain groups. Here H/L
forms one partner and A the other. dSASA reports two-sided buried area in Å²;
the CLI also reports its polarity components. SC is dimensionless.

```python
total, polar, hydrophobic, unclassified = arpeggia.dsasa_components(
    "complex.pdb", groups="H,L/A", probe_radius=1.4, n_points=100,
)
print(total, polar, hydrophobic, unclassified)
print(arpeggia.dsasa("complex.pdb", groups="H,L/A"))  # Total area only.
print(arpeggia.sc("complex.pdb", groups="H,L/A"))
```

```bash
arpeggia dsasa -i complex.pdb --groups "H,L/A" --probe-radius 1.4 --num-points 100
arpeggia sc -i complex.pdb --groups "H,L/A"
```

## Observed and declared sequences

`seq` uses residues with coordinates; `seqres` reads declared sequences, which
may include structurally unresolved residues.

```python
print(arpeggia.seq("structure.cif", model_num=0))
print(arpeggia.seqres("structure.cif"))
```

```bash
arpeggia seq structure.cif --model 0
arpeggia seqres structure.cif
```

For unaligned sequence strings, see the
[pairwise-alignment examples](https://github.com/y1zhou/arpeggia/blob/master/docs/sequence-alignment.md)
and [antibody-numbering examples](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md).

## Fit one region and evaluate another

Selections use author residue numbers. Here chain A establishes the fit and
chains B/C are evaluated after applying it. Each selection defaults independently
to all eligible residues, so provide both to fit and evaluate the same subset.

```python
result = arpeggia.rmsd(
    "reference.cif", "query.cif",
    superpose_residues="A:1-100,A:110-120", rmsd_residues="B,C", atoms="ca",
)
print(result.rmsd, result.core_rmsd)
```

```bash
arpeggia rmsd reference.cif query.cif --atoms ca \
  --superpose-residues "A:1-100,A:110-120" --rmsd-residues "B,C"
```

The [structure-comparison guide](https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md)
covers insertion codes, negative residue numbers and atom correspondence.
Enable `align_seqs=True` / `--align-seqs` to infer sequence correspondence;
then both selections use reference numbering. Use `refine_cycles=3` /
`--refine-cycles 3` for three rejection/refitting cycles (default: zero).

## Reuse a pairwise RMSD table for clustering

Collection calculations require exact atom correspondence and at least three
structures for clustering. Compute the pair table once, then reuse it for
different cluster counts:

```python
pairs = arpeggia.pairwise_rmsd("structures/", atoms="ca", num_threads=8)
pairs.write_parquet("pairwise_rmsd.parquet")
clusters = arpeggia.cluster_structs(pairwise_rmsd=pairs, num_clusters=3)
clusters.write_csv("clusters.csv")
```

```bash
arpeggia cluster-structs -i structures/ -o results/ --num-clusters 3 \
  --pairwise-rmsd --atoms ca --num-threads 8 -t parquet
```

The CLI reuses its existing pairwise table after checking schema and IDs; it
does not compare source contents or selections. Remove that table when inputs
or calculation settings change. See the
[cache contract](https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md#clustering).
