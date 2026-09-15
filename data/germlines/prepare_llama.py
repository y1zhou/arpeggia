"""Extract functional llama references from the pinned IMGT protein displays."""

import argparse
import hashlib
from html.parser import HTMLParser
from pathlib import Path


class ProteinTable(HTMLParser):
    """Read text cells, including styled sequence fragments, from IMGT tables."""

    def __init__(self) -> None:
        super().__init__()
        self.rows: list[list[str]] = []
        self.in_cell = False

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag == "tr":
            self.rows.append([])
        elif tag == "td":
            self.rows[-1].append("")
            self.in_cell = True

    def handle_endtag(self, tag: str) -> None:
        if tag == "td":
            self.in_cell = False

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.rows[-1][-1] += data


def main() -> None:
    """Retain functional gene/allele names, accessions and displayed IMGT gaps."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("v_display", type=Path)
    parser.add_argument("j_display", type=Path)
    args = parser.parse_args()
    sources = [
        (args.v_display, "e36c191249e7b67ea31da66e848e07ba09605a388c617842514f378d7cd005b7"),
        (args.j_display, "ec579819eef533d4031c857de310582d1663982b7d1ab23f9598273edd035e49"),
    ]
    for path, expected in sources:
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != expected:
            parser.error(f"{path.name} does not match the 2026-09-15 IMGT display")
        table = ProteinTable()
        table.feed(data.decode())
        for cells in table.rows:
            if len(cells) != 8 or cells[6] != "F":
                continue
            # These pinned displays have complete V frameworks; spaces separate
            # structural blocks or align J junctions, while dots are IMGT gaps.
            sequence = "".join(cells[7].split())
            print(f">{cells[4]}|{cells[3]}|Lama glama|F|{cells[5]}|protein-display")
            print(sequence)


if __name__ == "__main__":
    main()
