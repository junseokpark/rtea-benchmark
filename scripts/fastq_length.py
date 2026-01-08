#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
fastq_lenght.py

Compute per-read sequence lengths from FASTQ (.fastq/.fq) and gzipped FASTQ (.gz).
Output: <prefix>_lenghts.txt  (keeps original spelling for compatibility)
"""

import gzip
import os
import sys
from typing import TextIO


SYNTAX = """
------------------------------------------------------------------------------------
Usage: python fastq_lenght.py file.fastq[.gz]
result: .txt file same name as input name plus "_lenghts.txt"
------------------------------------------------------------------------------------
""".strip()


def open_maybe_gzip(path: str) -> TextIO:
    """Open plain text or .gz FASTQ as text."""
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "rt", encoding="utf-8", errors="replace", newline="")


def output_prefix(input_path: str) -> str:
    """Strip .gz (if present) then strip one more extension to create a stable prefix."""
    p = input_path
    if p.lower().endswith(".gz"):
        p = os.path.splitext(p)[0]  # remove .gz
    p = os.path.splitext(p)[0]      # remove .fastq/.fq/.txt/etc
    return p


def fastq_lengths(in_path: str) -> str:
    prefix = output_prefix(in_path)
    out_path = f"{prefix}_lenghts.txt"

    with open_maybe_gzip(in_path) as fin, open(out_path, "wt", encoding="utf-8", newline="") as fout:
        record_num = 0
        while True:
            header = fin.readline()
            if not header:
                break  # EOF

            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()

            if not (seq and plus and qual):
                raise ValueError(f"Truncated FASTQ record near record #{record_num + 1} in: {in_path}")

            header = header.rstrip("\n\r")
            seq = seq.rstrip("\n\r")
            plus = plus.rstrip("\n\r")
            qual = qual.rstrip("\n\r")

            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header at record #{record_num + 1}: {header!r}")
            if not plus.startswith("+"):
                raise ValueError(f"Invalid FASTQ '+' line at record #{record_num + 1}: {plus!r}")

            # We intentionally use sequence length (not quality).
            fout.write(f"{header}\t{len(seq)}\n")
            record_num += 1

    return out_path


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(SYNTAX)
        return 2

    in_path = argv[1]
    try:
        out_path = fastq_lengths(in_path)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"\n\tFile: {out_path} has been created...")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
