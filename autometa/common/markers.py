#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

Autometa Marker class consisting of various methods to annotate sequences with
marker sets depending on sequence set taxonomy
"""

import logging
import os

import pandas as pd

from autometa.common.external import hmmer


BASE_DIR = os.path.dirname(os.path.dirname(__file__))
MARKERS_DIR = os.path.join(BASE_DIR, "databases", "markers")

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
        return grouped_df.value_counts().unstack()
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
    dbdir: str,
    scans: str = None,
    out: str = None,
    force: bool = False,
    seed: int = 42,
    format: str = "wide",
    **hmmer_args: dict,
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
        Directory should contain hmmpressed marker genes database files.
    scans: str, optional
        Path to existing hmmscan table to filter by cutoffs
    out: str, optional
        Path to write annotated markers table.
    force : bool, optional
        Whether to overwrite existing `out` file path, by default False.
    format : str, optional
        * wide - returns wide dataframe of contig PFAM counts (default)
        * long - returns long dataframe of contig PFAM counts
        * list - returns list of pfams for each contig
        * counts - returns count of pfams for each contig
    hmmer_args: dict, optional
        additional keyword arguments to pass to `hmmer.hmmscan`

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
    orfs = os.path.realpath(orfs)
    kingdom = kingdom.lower()
    hmmdb = os.path.join(dbdir, f"{kingdom}.single_copy.hmm")
    cutoffs = os.path.join(dbdir, f"{kingdom}.single_copy.cutoffs")
    hmmscan_fname = ".".join([kingdom, "hmmscan.tsv"])
    scans = os.path.join(os.path.dirname(orfs), hmmscan_fname) if not scans else scans
    markers_fname = ".".join([kingdom, "markers.tsv"])
    out = os.path.join(os.path.dirname(orfs), markers_fname) if not out else out
    kingdom = kingdom.lower()

    if not os.path.exists(scans) or not os.path.getsize(scans):
        scans = hmmer.hmmscan(
            orfs=orfs, hmmdb=hmmdb, outfpath=scans, seed=seed, **hmmer_args,
        )

    if not os.path.exists(out) or not os.path.getsize(out):
        out = hmmer.filter_markers(
            infpath=scans, outfpath=out, cutoffs=cutoffs, orfs=orfs,
        )
    return load(fpath=out, format=format)


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
        "--force",
        help="Whether to overwrite existing provided annotations.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--seed", help=f"Seed to set random state for hmmscan.", default=42, type=int,
    )
    args = parser.parse_args()

    get(
        kingdom=args.kingdom,
        orfs=args.orfs,
        dbdir=args.dbdir,
        scans=args.hmmscan,
        out=args.out,
        force=args.force,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
