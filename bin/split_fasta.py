#!/usr/bin/env python

# Modified by Chase Clark 2023
# Copyright 2006-2017,2020 by Peter Cock.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
#
# This module is for reading and writing FASTA format files as SeqRecord
# objects.  The code is partly inspired  by earlier Biopython modules,
# Bio.Fasta.* and the now removed module Bio.SeqIO.FASTA

# python dependencies
import bz2
import gzip
from contextlib import contextmanager
from enum import Enum, auto
from pathlib import Path
from typing import TextIO
import math
import argparse

parser = argparse.ArgumentParser(description="Split fasta file into x-files")

parser.add_argument(
    "--input",
    metavar="filepath",
    help="fasta file path",
    required=True,
)
parser.add_argument(
    "--outdir",
    metavar="filepath",
    help="path split files will be written to",
    required=True,
)
parser.add_argument(
    "--splits",
    metavar="int",
    help="number of files to split input into",
    required=True,
)


def SimpleFastaParser(handle):
    # Skip any text before the first record (e.g. blank lines, comments)
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break
    else:
        # no break encountered - probably an empty file
        return
    # Main logic
    # Note, remove trailing whitespace, and any internal spaces
    # (and any embedded \r which are possible in mangled files
    # when not opened in universal read lines mode)
    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())
    yield title, "".join(lines).replace(" ", "").replace("\r", "")


class Compression(Enum):
    bzip2 = auto()
    gzip = auto()
    xz = auto()
    uncompressed = auto()


def is_compressed(filepath: Path) -> Compression:
    with open(filepath, "rb") as f:
        signature = f.peek(8)[:8]
        if tuple(signature[:2]) == (0x1F, 0x8B):
            return Compression.gzip
        elif tuple(signature[:3]) == (0x42, 0x5A, 0x68):
            return Compression.bzip2
        else:
            return Compression.uncompressed


@contextmanager
def open_file(filepath: Path) -> TextIO:
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        f = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        f = bz2.open(filepath, "rt")
    else:
        f = open(filepath, "r")
    try:
        yield f
    finally:
        f.close()


def main():

    args = parser.parse_args()

    input = Path(args.input)
    if not input.exists():
        raise FileExistsError

    outdir = Path(args.outdir)
    if not outdir.exists():
        raise FileExistsError

    limit = int(args.splits)
    records = 0
    with open_file(input) as handle:
        for record in SimpleFastaParser(handle):
            records += 1

    records_per_file = math.ceil(records / limit)
    outfiles = (
        gzip.open(Path(outdir, f"{input.stem}_{i}.gz"), "wb", compresslevel=6)
        for i in range(0, limit)
    )
    outfile = outfiles.__next__()
    file_count = 0
    n = 0
    with open_file(input) as handle:
        for line in SimpleFastaParser(handle):
            if line:
                if n == records_per_file and file_count < limit:
                    n = 0
                    outfile.close()
                    outfile = outfiles.__next__()
                    print(f"file {file_count} of {limit}", end="\r")
                outfile.write(str.encode(f">{line[0]}\n{line[1]}\n"))
                n += 1
        outfile.close()


if __name__ == "__main__":
    main()
