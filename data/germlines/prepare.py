"""Reproduce the attributed IMGT protein subset from the pinned source download."""

import argparse
import hashlib
import sys
from pathlib import Path


def main() -> None:
    """Filter functional antibody V/J records without discarding source metadata."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    data = args.source.read_bytes()
    expected = "3cb6b0b8cb8940b3b2a9b105771a6a74aa67c06e3ca39eaea0d2030c90e7efd0"
    if hashlib.sha256(data).hexdigest() != expected:
        parser.error("source does not match IMGT release 202636-7")
    for record in data.decode().split(">")[1:]:
        header, *lines = record.splitlines()
        fields = header.split("|")
        sequence = "".join(lines)
        if (
            fields[1].startswith(("IGHV", "IGKV", "IGLV", "IGHJ", "IGKJ", "IGLJ"))
            and fields[2].startswith(
                (
                    "Homo sapiens",
                    "Mus musculus",
                    "Vicugna pacos",
                    "Rattus norvegicus",
                    "Oryctolagus cuniculus",
                )
            )
            and fields[3].strip("()[]") == "F"
            and "*" not in sequence
        ):
            sys.stdout.write(f">{header}\n{sequence}\n")


if __name__ == "__main__":
    main()
