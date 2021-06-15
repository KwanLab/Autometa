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
import logging
import pandas as pd

logger = logging.getLogger(__name__)


def download(dataset: str, file: str, out: str) -> None:

    """Downloads the files specified in a dictionary.

    Parameters
    ----------
    dataset : str
        specifies the dataset that the user would like to download from, if any
    file : str
        specifies the file that the user would like to download, if any
    output : str
        directory path where the user wants to download the file(s)

    Returns
    -------
    None
        download is completed through gdown

    """
    # create the dataframe here
    index = pd.read_csv("gdown_fileIDs.csv", dtype=str)

    # retrieve only the community that the user wants (user must specify a community)
    target_dataset = index.query(f'dataset == "{dataset}"')

    target_files = target_dataset.query(f'file == "{file}"')

    targets = target_files.to_dict()

    key_list = [*targets["file"]]
    for key in key_list:
        # retrieve file ids from targets
        file_id = targets["file_id"][key]

        # construct file id into a url to put into gdown
        url = f"https://drive.google.com/uc?id={file_id}"

        # download the specified file with gdown
        gdown.download(url, out)


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Download a simulated community file from google drive to a specified directory"
    )
    parser.add_argument(
        "--community",
        help="specify a simulated community size to download from",
        choices=[
            "78.125Mbp",
            "156.25Mbp",
            "312.5Mbp",
            "625Mbp",
            "1250Mbp",
            "2500Mbp",
            "5000Mbp",
            "10000Mbp",
        ],
        required=True,
    )
    parser.add_argument(
        "--file",
        help="specify a file to download",
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
        ],
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "--output",
        help="specify the full filepath for your downloaded file (including filename)",
        required=True,
        nargs="+",
    )
    args = parser.parse_args()

    # if I add in "all" functionality, add it here

    community = args.community
    filetypes = args.file
    filepaths = args.output

    if len(filetypes) != len(filepaths):
        logger.warning(
            "The number of files specified and the number of output paths specified do not match. The program will use the shorter list. You might not get all the files you wanted."
        )
    for filetype, filepath in zip(filetypes, filepaths):
        if (
            os.path.splitext(filetype)[-1].lower()
            != os.path.splitext(filepath)[-1].lower()
        ):
            logger.warning(
                "The file extension doesn't match on the file and the output. The file will still be saved under the user-specified output, but it might be difficult to open."
            )
        download(dataset=community, file=filetype, out=filepath)


if __name__ == "__main__":
    main()
