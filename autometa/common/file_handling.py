# Copied from SocialGene
# The MIT License (MIT)
#
# Copyright (c) 2021 Chase M Clark
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# python dependencies
import bz2
import gzip
import lzma
from contextlib import contextmanager
from enum import Enum, auto
from pathlib import Path
from typing import TextIO
import tarfile
import shutil
import logging as log


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
        elif tuple(signature[:7]) == (0xFD, 0x37, 0x7A, 0x58, 0x5A, 0x00, 0x00):
            return Compression.xz
        else:
            return Compression.uncompressed


def gunzip(filepath: Path) -> None:
    if is_compressed(filepath).name == "gzip":
        # remove ".gz" for the new filepath
        new_hmm_path = Path(filepath.parents[0], filepath.stem)
        log.info(f"Start decompressing: {str(Path(filepath).stem)}")
        # open gz, decompress, write back out to new file
        with gzip.open(filepath, "rb") as f_in:
            with open(new_hmm_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        log.info(f"Finish decompressing: {str(Path(filepath).stem)}")


@contextmanager
def open_file(filepath: Path) -> TextIO:
    filepath_compression = is_compressed(filepath)
    if filepath_compression == Compression.gzip:
        f = gzip.open(filepath, "rt")
    elif filepath_compression == Compression.bzip2:
        f = bz2.open(filepath, "rt")
    elif filepath_compression == Compression.xz:
        f = lzma.open(filepath, "rt")
    else:
        f = open(filepath, "r")
    try:
        yield f
    finally:
        f.close()


def check_if_tar(filepath):
    return tarfile.is_tarfile(filepath)
