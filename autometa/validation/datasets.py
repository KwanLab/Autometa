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
    community_type: str, community_size: list, file_names: list, dir_path: str
) -> None:

    """Downloads the files specified in a dictionary.

    Parameters
    ----------
    community_type : str
        specifies the type of dataset that the user would like to download from
    community_size : list
        specifies the size of dataset that the user would like to download
    file_names : list
        specifies the file(s) that the user would like to download
    dir_path : str
        dir_path_size path where the user wants to download the file(s)

    Returns
    -------
    None
        download is completed through gdown

    """

    if community_type == "synthetic" or community_type == "all":
        raise NotImplementedError

    for each_community_size in community_size:
        df = pd.read_csv("gdown_fileIDs.csv", dtype=str)
        df = df.query(f'dataset == "{each_community_size}"')
        dir_path_size = os.path.join(dir_path, each_community_size)
        # make a new directory
        if not os.path.exists(dir_path_size):
            os.mkdir(dir_path_size)

        for file_name in file_names:
            file_id = df.query(f'file == "{file_name}"')["file_id"].to_list()[0]
            dir_path_final = os.path.join(dir_path_size, file_name)
            url = f"https://drive.google.com/uc?id={file_id}"

            gdown.download(url, dir_path_final)


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Download a simulated community file from google drive to a specified dir_path_size"
    )
    parser.add_argument(
        "--community_type",
        help="specify synthetic or simulated communities (currently only simulated is available)",
        choices=[
            "synthetic",
            "simulated",
            "all",
        ],
        required=True,
    )
    parser.add_argument(
        "--community_size",
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
        "--file_names",
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
        "--dir_path",
        help="specify a folder to start the download (several directories will be generated within this folder)",
        required=True,
    )
    args = parser.parse_args()

    community_type = args.community_type
    if "all" in args.community_size:
        community_size = (
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
        community_size = args.community_size
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
    dir_path = args.dir_path

    if not internet_is_connected():
        logger.error(
            "No internet connection detected. Please confirm connection. Downloader will still attempt to run. (Internet check may incorrectly fail where google.com is not reachable/ping-able (e.g. China))"
        )

    if not os.path.exists(dir_path):
        logger.error(
            "The dir_path you specified could not be found. Please select a folder that exists."
        )
        sys.exit("Specified path does not exist!")

    download(
        community_type=community_type,
        community_size=community_size,
        file_names=file_names,
        dir_path=dir_path,
    )


if __name__ == "__main__":
    main()
