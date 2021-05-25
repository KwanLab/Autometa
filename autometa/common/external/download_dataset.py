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

pulling data from google drive folder with simulated or synthetic communities
"""


import gdown
import os
import logging

logger = logging.getLogger(__name__)


def download_dataset(dataset, out_dirpath):
    # provide list of database options as a dictionary with file_ids from google
    simulated = {
        "test": "1fy3M7RnS_HGSQVKidCy-rAwXuxldyOOv",
        "78": "15CB8rmQaHTGy7gWtZedfBJkrwr51bb2y",
        "156": "13bkwFBIUhdWVWlAmVCimDODWF-7tRxgI",
        "312": "1qyAu-m6NCNuVlDFFC10waOD28j15yfV-",
        "625": "1FgMXSD50ggu0UJbZd1PM_AvLt-E7gJix",
        "1250": "1KoxwxBAYcz8Xz9H2v17N9CHOZ-WXWS5m",
        "2500": "1wKZytjC4zjTuhHdNUyAT6wVbuDDIwk2m",
        "5000": "1IX6vLfBptPxhL44dLa6jePs-GRw2XJ3S",
        "10000": "1ON2vxEWC5FHyyPqlfZ0znMgnQ1fTirqG",
    }

    # construct file id into a url to put into gdown
    file_id = simulated[dataset]
    url = f"https://drive.google.com/uc?id={file_id}"
    filename = f"{dataset}_metagenome.fna.gz"
    out_fpath = os.path.join(out_dirpath, filename)

    # download the specified file with gdown
    gdown.download(url, out_fpath)


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        prog="autometa-download-dataset",
        description="Download a simulated community file from google drive to a specified directory",
    )
    parser.add_argument(
        "--dataset",
        help="specify a size of simulated community in megabase pairs",
        choices=["78", "156", "312", "625", "1250", "2500", "5000", "10000", "test"],
        required=True,
    )
    parser.add_argument(
        "--out_dirpath",
        help="specify the directory to download the file",
        required=True,
    )
    args = parser.parse_args()

    download_dataset(args.dataset, args.out_dirpath)


if __name__ == "__main__":
    main()
