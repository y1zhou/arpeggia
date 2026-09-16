"""Exercise an installed wheel without importing the source tree."""

from pathlib import Path

import arpeggia

sequences = arpeggia.seq(str(Path(__file__).parents[2] / "test-data" / "1ubq.pdb"))
expected = [
    (
        "A",
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    )
]
if sequences != expected:
    raise SystemExit(f"unexpected observed sequence: {sequences!r}")

# References and metadata must work from an installed wheel without source data.
antibody = arpeggia.number_antibody(
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGR"
    "FTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRGGYFDYWGQGTLVTVSS"
)
if antibody.chain != "H" or antibody.v_match is None or antibody.j_match is None:
    raise SystemExit(
        "installed wheel failed antibody numbering or offline V/J matching"
    )
alignment = arpeggia.align_antibodies([antibody])
if alignment.aligned_sequences != [antibody.sequence]:
    raise SystemExit("installed wheel failed numbered antibody alignment")
