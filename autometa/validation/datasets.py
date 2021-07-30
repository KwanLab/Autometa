#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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

pulling data from google drive dataset with simulated or synthetic communities
"""


import gdown
import os
import sys
import logging
import pandas as pd

from autometa.common.utilities import internet_is_connected

logger = logging.getLogger(__name__)


def download(
    community_type: str, community_sizes: list, file_names: list, dir_path: str
) -> None:

    """Downloads the files specified in a dictionary.

    Parameters
    ----------
    community_type : str
        specifies the type of dataset that the user would like to download from
    community_sizes : list
        specifies the size of dataset that the user would like to download
    file_names : list
        specifies the file(s) that the user would like to download
    dir_path : str
        output path where the user wants to download the file(s)

    Returns
    -------
    None
        download is completed through gdown

    """

    if community_type == "synthetic" or community_type == "all":
        raise NotImplementedError

    # points to csv file on google drive
    df = pd.read_csv(
        "https://drive.google.com/uc?id=148fUO7jocoNOBUl2K4bCfjsbd42QxCzX",
        dtype=str,
        index_col=["dataset", "file"],
    )

    for community_size in community_sizes:
        community_size_outdir = os.path.join(dir_path, community_size)
        # make a new directory
        if not os.path.exists(community_size_outdir):
            os.makedirs(community_size_outdir)

        for file_name in file_names:
            file_id = df.loc[(community_size, file_name), "file_id"]
            file_id_filepath = os.path.join(community_size_outdir, file_name)
            url = f"https://drive.google.com/uc?id={file_id}"

            gdown.download(url, file_id_filepath)


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Download a simulated community file from google drive to a specified output directory"
    )
    parser.add_argument(
        "--community-type",
        help="specify synthetic or simulated communities (currently only simulated is available)",
        choices=[
            "synthetic",
            "simulated",
            "all",
        ],
        required=True,
    )
    parser.add_argument(
        "--community-sizes",
        help="specify a community size to download from",
        choices=[
            "78Mbp",
            "156Mbp",
            "312Mbp",
            "625Mbp",
            "1250Mbp",
            "2500Mbp",
            "5000Mbp",
            "10000Mbp",
            "all",
        ],
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "--file-names",
        help="specify a file name to download",
        choices=[
            "README.md",
            "reference_assignments.tsv.gz",
            "metagenome.fna.gz",
            "master.tsv.gz",
            "control_reads.tsv.gz",
            "control_contigs.tsv.gz",
            "unclustered_recruitment.tsv.gz",
            "binning.tsv.gz",
            "taxonomy.tsv.gz",
            "lengths.tsv.gz",
            "coverages.tsv.gz",
            "gc_content.tsv.gz",
            "kmers.embedded.tsv.gz",
            "kmers.tsv.gz",
            "markers.tsv.gz",
            "Bacteria.fna.gz",
            "orfs.faa.gz",
            "metagenome.filtered.fna.gz",
            "hmmscan.tsv.gz",
            "forward_reads.fastq.gz",
            "reverse_reads.fastq.gz",
            "all",
        ],
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "--dir-path",
        help="specify a folder to start the download (several directories will be generated within this folder)",
        required=True,
    )
    parser.add_argument(
        "--host",
        help="IP address to ping when checking internet connectivity. Note: Will attempt to connect to port 53 on host address (Default is google.com)",
        default="8.8.8.8",
    )
    args = parser.parse_args()

    if "all" in args.community_sizes:
        community_sizes = (
            "78Mbp",
            "156Mbp",
            "312Mbp",
            "625Mbp",
            "1250Mbp",
            "2500Mbp",
            "5000Mbp",
            "10000Mbp",
        )
    else:
        community_sizes = args.community_sizes
    if "all" in args.file_names:
        file_names = (
            "README.md",
            "reference_assignments.tsv.gz",
            "metagenome.fna.gz",
            "master.tsv.gz",
            "control_reads.tsv.gz",
            "control_contigs.tsv.gz",
            "unclustered_recruitment.tsv.gz",
            "binning.tsv.gz",
            "taxonomy.tsv.gz",
            "lengths.tsv.gz",
            "coverages.tsv.gz",
            "gc_content.tsv.gz",
            "kmers.embedded.tsv.gz",
            "kmers.tsv.gz",
            "markers.tsv.gz",
            "Bacteria.fna.gz",
            "orfs.faa.gz",
            "metagenome.filtered.fna.gz",
            "hmmscan.tsv.gz",
            "forward_reads.fastq.gz",
            "reverse_reads.fastq.gz",
        )
    else:
        file_names = args.file_names

    if not internet_is_connected(host=args.host):
        logger.error(
            "No internet connection detected (couldn't ping google.com at IP 8.8.8.8). Please confirm connection. Downloader will still attempt to run. (Ping a custom IP address with --host argument)"
        )

    download(
        community_type=args.community_type,
        community_sizes=community_sizes,
        file_names=file_names,
        dir_path=args.dir_path,
    )


if __name__ == "__main__":
    main()
