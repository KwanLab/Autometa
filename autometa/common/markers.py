#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Autometa Marker class consisting of various methods to annotate sequences with
marker sets depending on sequence set taxonomy
"""

import logging
import os
import sys

import pandas as pd

from autometa.common.external import hmmscan
from autometa.config.utilities import DEFAULT_CONFIG


MARKERS_DIR = DEFAULT_CONFIG.get("databases", "markers")
# For cases where autometa has not been configured, attempt to find the markers via source
MARKERS_DIR = (
    MARKERS_DIR if not "None" in MARKERS_DIR else MARKERS_DIR.replace("None", ".")
)

logger = logging.getLogger(__name__)


def load(fpath, format="wide"):
    """Read markers table into specified `format`.

    Parameters
    ----------
    fpath : str
        </path/to/`kingdom`.markers.tsv>
    format : str, optional
        * wide - index=contig, cols=[domain sacc,..] (default)
        * long - index=contig, cols=['sacc','count']
        * list - {contig:[sacc,...],...}
        * counts - {contig:len([sacc,...]), ...}

    Returns
    -------
    pd.DataFrame or dict
        * wide - index=contig, cols=[domain sacc,..] (default)
        * long - index=contig, cols=['sacc','count']
        * list - {contig:[sacc,...],...}
        * counts - {contig:len([sacc,...]), ...}

    Raises
    -------
    FileNotFoundError
        Provided `fpath` does not exist
    ValueError
        Provided `format` is not in choices:
        choices = wide, long, list or counts

    """
    if not os.path.exists(fpath):
        raise FileNotFoundError(fpath)
    df = pd.read_csv(fpath, sep="\t", index_col="contig")
    grouped_df = df.groupby("contig")["sacc"]
    if format == "wide":
        return grouped_df.value_counts().unstack().convert_dtypes()
    elif format == "long":
        return grouped_df.value_counts().reset_index(level=1, name="count")
    elif format == "list":
        return {contig: markers.tolist() for contig, markers in list(grouped_df)}
    elif format == "counts":
        return grouped_df.count().to_dict()
    else:
        params = ["wide", "long", "list", "counts"]
        err_msg = f"{format} is not a supported format.\n\tSupported formats: {params}"
        # TODO: Write Marker specific AutometaException
        raise ValueError(err_msg)


def get(
    kingdom: str,
    orfs: str,
    hmmdb: str = None,
    cutoffs: str = None,
    dbdir: str = MARKERS_DIR,
    scans: str = None,
    out: str = None,
    force: bool = False,
    cpus: int = 8,
    parallel: bool = True,
    gnu_parallel: bool = False,
    seed: int = 42,
) -> pd.DataFrame:
    """Retrieve contigs' markers from markers database that pass cutoffs filter.

    Parameters
    ----------
    kingdom: str
        kingdom to annotate markers
        choices = ['bacteria', 'archaea']
    orfs: str
        Path to amino-acid ORFs file
    dbdir:
        Optional directory containing hmmdb and cutoffs files
    hmmdb:
        Path to marker genes database file, previously hmmpressed.
    cutoffs:
        Path to marker genes cutoff tsv.
    scans: str, optional
        Path to existing hmmscan table to filter by cutoffs
    out: str, optional
        Path to write annotated markers table.
    force : bool, optional
        Whether to overwrite existing `out` file path, by default False.
    cpus: int, optional
        Number of cores to use if running in parallel, by default all available.
    parallel : bool, optional
        Whether to run hmmscan using its parallel option, by default True.
    gnu_parallel : bool, optional
        Whether to run hmmscan using gnu parallel, by default False.
    seed: int, optional
        Seed to pass into hmmscan for determinism, by default 42.

    Returns
    -------
    pd.Dataframe or dict
        * wide - pd.DataFrame(index_col=contig, columns=[PFAM,...])
        * long - pd.DataFrame(index_col=contig, columns=['sacc','count'])
        * list - {contig:[pfam,pfam,...],contig:[...],...}
        * counts - {contig:count, contig:count,...}

    Raises
    -------
    ValueError
        Why the exception is raised.
    """
    kingdom = kingdom.lower()
    # if dbdir == MARKERS_DIR and hmmdb/cutoffs not set, use single dbdir
    # else dbdir was set, so use it for both hmmdb and cutoffs
    if dbdir == MARKERS_DIR:
        if hmmdb is None:
            hmmdb = os.path.join(dbdir, f"{kingdom}.single_copy.hmm")
        if cutoffs is None:
            cutoffs = os.path.join(dbdir, f"{kingdom}.single_copy.cutoffs")
    else:
        hmmdb = os.path.join(dbdir, f"{kingdom}.single_copy.hmm")
        cutoffs = os.path.join(dbdir, f"{kingdom}.single_copy.cutoffs")
    hmmscan_fname = f"{kingdom}.hmmscan.tsv"
    scans = (
        os.path.join(os.path.dirname(os.path.abspath((orfs))), hmmscan_fname)
        if not scans
        else scans
    )
    markers_fname = f"{kingdom}.markers.tsv"
    out = (
        os.path.join(os.path.dirname(os.path.abspath((orfs))), markers_fname)
        if not out
        else out
    )
    kingdom = kingdom.lower()

    if not os.path.exists(scans) or not os.path.getsize(scans):
        scans = hmmscan.run(
            orfs=orfs,
            hmmdb=hmmdb,
            outfpath=scans,
            cpus=cpus,
            force=force,
            parallel=parallel,
            gnu_parallel=gnu_parallel,
            seed=seed,
        )
        logger.debug(f"{kingdom} marker scan finished")
    if not os.path.exists(out) or not os.path.getsize(out):
        logger.debug(f"filtering {scans} for marker genes passing cutoffs")
        df = hmmscan.filter_tblout_markers(
            infpath=scans,
            outfpath=out,
            cutoffs=cutoffs,
            orfs=orfs,
            force=force,
        )
        logger.debug(f"marker gene cutoff filter finished: {out}")
        grouped_df = df.groupby("contig")["sacc"]
        return grouped_df.value_counts().unstack().convert_dtypes()
    else:
        return load(fpath=out, format="wide")


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Annotate ORFs with kingdom-specific marker information",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--orfs",
        help="Path to a fasta file containing amino acid sequences of open reading frames",
    )
    parser.add_argument(
        "--kingdom",
        help="kingdom to search for markers",
        choices=["bacteria", "archaea"],
        default="bacteria",
    )
    parser.add_argument(
        "--hmmscan",
        help="Path to hmmscan output table containing the respective `kingdom` single-copy marker annotations.",
    )
    parser.add_argument(
        "--out",
        help="Path to write filtered annotated markers corresponding to `kingdom`.",
    )
    parser.add_argument(
        "--dbdir",
        help="Path to directory containing the single-copy marker HMM databases.",
        default=MARKERS_DIR,
    )
    parser.add_argument(
        "--hmmdb",
        help="Path to single-copy marker HMM databases.",
    )
    parser.add_argument(
        "--cutoffs",
        help="Path to single-copy marker cutoff tsv.",
    )
    parser.add_argument(
        "--force",
        help="Whether to overwrite existing provided annotations.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--parallel",
        help="Whether to use hmmscan parallel option.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--gnu-parallel",
        help="Whether to run hmmscan using GNU parallel.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--cpus",
        help=f"Number of cores to use for parallel execution.",
        default=8,
        type=int,
    )
    parser.add_argument(
        "--seed",
        help="Seed to set random state for hmmscan.",
        default=42,
        type=int,
    )
    args = parser.parse_args()

    try:
        get(
            kingdom=args.kingdom,
            orfs=args.orfs,
            hmmdb=args.hmmdb,
            dbdir=args.dbdir,
            cutoffs=args.cutoffs,
            scans=args.hmmscan,
            out=args.out,
            force=args.force,
            cpus=args.cpus,
            parallel=args.parallel,
            gnu_parallel=args.gnu_parallel,
            seed=args.seed,
        )
    except AssertionError as err:
        logger.warn(err)
        sys.exit(204)


if __name__ == "__main__":
    main()
