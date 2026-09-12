# Python recipes

Install Arpeggia as described in the
[README](https://github.com/y1zhou/arpeggia/blob/master/README.md#installation).
These examples use your own `structure.pdb`; interface examples require the
specified chains. All tabular results are Polars DataFrames.

## Contacts and hydrogen bonds

```python
import arpeggia
import polars as pl

contacts = arpeggia.contacts("structure.pdb")
print(contacts.group_by("interaction").len())
hydrogen_bonds = contacts.filter(
    pl.col("interaction").is_in(["HydrogenBond", "WeakHydrogenBond"])
)
residue_pairs = hydrogen_bonds.group_by([
    "model", "from_chain", "from_resi", "from_insertion", "from_resn",
    "to_chain", "to_resi", "to_insertion", "to_resn",
]).len()
print(residue_pairs)

contacts.write_csv("contacts.csv")
contacts.write_parquet("contacts.parquet")
```

Residue identity includes model, chain, residue number and insertion code;
residue numbers alone can merge distinct residues. Counts above are contact
rows, not unique atom pairs: a pair can have multiple interaction types.
Hydrogen-bond results depend on hydrogen evidence and the selected protonation
policy; see the
[scientific conventions](https://github.com/y1zhou/arpeggia/blob/master/docs/scientific-conventions.md).

## Interface residues

```python
contacts = arpeggia.contacts("structure.pdb", groups="A/B")
identity = ["chain", "resi", "insertion", "resn"]
interface_residues = pl.concat([
    contacts.select("model", *[
        pl.col(f"{side}_{field}").alias(field) for field in identity
    ])
    for side in ("from", "to")
]).unique()
print(interface_residues)
```

Combine both endpoints: `from` and `to` can follow interaction roles rather
than chain-group order. Filter the resulting `chain` column for one partner. See the
[chain-group syntax](https://github.com/y1zhou/arpeggia/blob/master/README.md#chain-groups-specification)
for larger interfaces.

## Solvent exposure

```python
residue_sasa = arpeggia.sasa("structure.pdb", level="residue")
print(residue_sasa.sort("sasa").head(10))
relative_sasa = arpeggia.relative_sasa("structure.pdb")
print(relative_sasa)
```

Residue SASA sums atomic accessible areas in Å². `relative_sasa()` divides each
standard residue's area by its reference maximum; averaging atomic SASA would
measure a different quantity. Use `chains="A,B"` to restrict either calculation.

## SAP scores

```python
import arpeggia
import polars as pl

# Calculate SAP scores for an antibody
sap = arpeggia.sap_score("antibody.pdb", level="residue")

# Find aggregation-prone residues (high positive SAP)
aggregation_prone = sap.filter(pl.col("sap_score") > 0.5)
print(f"Found {len(aggregation_prone)} aggregation-prone residues")

# SAP for only heavy and light chains
sap_hl = arpeggia.sap_score("antibody.pdb", chains="H,L")
print(f"Analyzed {len(sap_hl)} residues in H and L chains")

# Sort by SAP score to find hotspots
hotspots = sap.sort("sap_score", descending=True).head(10)
print("Top 10 aggregation hotspots:")
print(hotspots)
```

## Shape complementarity

```python
import arpeggia

# Calculate Shape Complementarity between antibody and antigen
sc_score = arpeggia.sc("antibody_antigen.pdb", groups="H,L/A")
print(f"Shape Complementarity: {sc_score:.3f}")

# Typical SC values:
# - Good fit: 0.65-0.80
# - Average fit: 0.50-0.65
# - Poor fit: < 0.50

# SC between different chain groups
sc_hl = arpeggia.sc("antibody.pdb", groups="H/L")
print(f"VH-VL interface SC: {sc_hl:.3f}")
```

## Further usage

- [Sequence alignment and sequence-aware RMSD](https://github.com/y1zhou/arpeggia/blob/master/docs/sequence-alignment.md)
- [Structure selections, pairwise RMSD and clustering](https://github.com/y1zhou/arpeggia/blob/master/docs/structure-comparison.md)
- [Antibody numbering, germline matching and imputation](https://github.com/y1zhou/arpeggia/blob/master/docs/antibody-numbering.md)
- [Building and testing](https://github.com/y1zhou/arpeggia/blob/master/BUILD.md)
