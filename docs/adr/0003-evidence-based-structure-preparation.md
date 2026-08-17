# Use evidence-based structure preparation

Arpeggia reports definitive covalent contacts only when bonding information is
extractable from the input structure. Contacts meeting the existing
distance-based covalent criterion without that evidence are reported separately
as `PotentialCovalent`. Hydrogen-bond geometry uses only hydrogens bonded to the
specific donor; Arpeggia does not automatically protonate structures and warns
when possible donors lack hydrogens.

Histidine charge is controlled by `AllCharged`, `Heuristic`, or `ExplicitOnly`.
`AllCharged` is the default and preserves Arpeggio's unconditional
positive-ionisable typing. `Heuristic` uses explicit evidence first and otherwise
applies a documented pH-dependent estimate; estimated charge produces potential,
not definitive, interaction categories. `ExplicitOnly` accepts positive atom-site
formal charge, explicit doubly protonated ring hydrogens, or guarded recognized
protonation aliases and leaves all other histidines unknown. Each mode emits
structured warnings for assumptions or unknown protonation states.

The implementation carries searchable `[WARNING]` comments explaining that
intrinsic-pKa histidine heuristics cannot model the local molecular environment
and that selecting one alternate conformer does not represent an
occupancy-weighted conformational ensemble.

All calculations select one alternate conformer consistently: an unambiguous
blank/default conformer, otherwise the highest-occupancy conformer with `A` as a
deterministic tie-breaker. Automatic selection emits a warning rather than
silently mixing mutually exclusive coordinates.

The existing `seq` operation returns coordinate-observed polymer sequences and
excludes solvent and ligands. The separate `seqres` API exposes declared PDB
`SEQRES` or mmCIF entity sequences, including unresolved polymer residues. A
narrow metadata layer reads PDB `SSBOND`, `LINK`, `CONECT`, and `SEQRES`, plus
mmCIF `_struct_conn`, `_entity_poly_seq`, and `_pdbx_poly_seq_scheme`, while
`pdbtbx` remains the coordinate parser.

`seqres` remains chain-oriented like `seq`: it returns one declared sequence per
chain in declaration order, repeats a shared mmCIF entity sequence for each
mapped chain, and includes declared chains without coordinates. Unrecognized
monomers become `X` with a structured warning.
