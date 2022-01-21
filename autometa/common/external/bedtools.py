#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
Script containing wrapper functions for bedtools.
"""


import logging
import os
import subprocess

import pandas as pd

from autometa.common.exceptions import TableFormatError


logger = logging.getLogger(__name__)


def genomecov(ibam: str, out: str, force: bool = False) -> str:
    """Run bedtools genomecov with input `ibam` and `lengths` to retrieve
    metagenome coverages.

    Parameters
    ----------
    ibam : str
        </path/to/indexed/BAM/file.ibam>. Note: BAM *must* be sorted by position.

    out : str
        </path/to/alignment.bed>
        The bedtools genomecov output is a tab-delimited file with the following columns:
        1. Chromosome
        2. Depth of coverage
        3. Number of bases on chromosome with that coverage
        4. Size of chromosome
        5. Fraction of bases on that chromosome with that coverage
        See also: http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html

    force : bool
        force overwrite of `out` if it already exists (default is False).

    Returns
    -------
    str
        </path/to/alignment.bed>

    Raises
    -------
    FileExistsError
        `out` file already exists and force is False
    OSError
        Why the exception is raised.

    """
    cmd = f"bedtools genomecov -ibam {ibam}"
    if os.path.exists(out) and not force:
        logger.debug(f"{out} already exists. skipping...")
        return out
    with open(os.devnull, "w") as stderr, open(out, "w") as stdout:
        retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
    if retcode or not os.path.exists(out) or not os.path.getsize(out):
        raise ChildProcessError(f"bedtools failed: {cmd}")
    return out


def parse(bed: str, out: str = None, force: bool = False) -> pd.DataFrame:
    """Calculate coverages from bed file.

    Parameters
    ----------
    bed : str
        </path/to/file.bed>

    out : str
        if provided will write to `out`. I.e. </path/to/coverage.tsv>

    force : bool
        force overwrite of `out` if it already exists (default is False).

    Returns
    -------
    pd.DataFrame
        index='contig', col='coverage'

    Raises
    -------
    ValueError
        `out` incorrectly formatted to be read as pandas DataFrame.
    FileNotFoundError
        `bed` does not exist

    """
    if out and os.path.exists(out) and os.path.getsize(out):
        try:
            cols = ["contig", "coverage"]
            return pd.read_csv(out, sep="\t", usecols=cols, index_col="contig")
        except ValueError:
            raise TableFormatError(out)
    if not os.path.exists(bed):
        raise FileNotFoundError(bed)
    names = ["contig", "depth", "bases", "length", "depth_fraction"]
    df = pd.read_csv(bed, sep="\t", names=names, index_col="contig")
    criterion1 = df.depth != 0
    criterion2 = df.index != "genome"
    df = df[criterion1 & criterion2]
    df = df.assign(depth_product=lambda x: x.depth * x.bases)
    dff = df.groupby("contig")["depth_product", "bases"].sum()
    dff = dff.assign(coverage=lambda x: x.depth_product / x.bases)
    if out and (not os.path.exists(out) or (os.path.exists(out) and force)):
        dff.to_csv(out, sep="\t", index=True, header=True)
        logger.debug(f"{out} written")
    msg = (
        f"{os.path.basename(out)} shape: {dff.shape}" if out else f"shape: {dff.shape}"
    )
    logger.debug(msg)
    return dff[["coverage"]]


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="%(asctime)s : %(name)s : %(levelname)s : %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    parser = argparse.ArgumentParser(
        description="Compute genome coverage from sorted BAM file"
    )
    parser.add_argument(
        "--ibam", metavar="filepath", help="Path to sorted alignment.bam", required=True
    )
    parser.add_argument(
        "--bed",
        metavar="filepath",
        help="Path to write alignment.bed; tab-delimited cols=[contig,length]",
        required=True,
    )
    parser.add_argument(
        "--output",
        metavar="filepath",
        help="Path to output coverage.tsv",
        required=True,
    )
    parser.add_argument(
        "--force-bed", help="force overwrite `bed`", action="store_true", default=False
    )
    parser.add_argument(
        "--force-cov",
        help="force overwrite `--output`",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()
    bed = genomecov(ibam=args.ibam, out=args.bed, force=args.force_bed)
    parse(bed=bed, out=args.output, force=args.force_cov)


if __name__ == "__main__":
    main()
