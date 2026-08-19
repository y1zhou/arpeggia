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
