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
Script containing wrapper functions for samtools
"""


import logging
import os
import shutil
import subprocess
import tempfile

import pandas as pd
import multiprocessing as mp

logger = logging.getLogger(__name__)


def sort(sam, bam, cpus=mp.cpu_count()):
    """
    Views then sorts sam file by leftmost coordinates and outputs to bam.

    Parameters
    ----------
    sam : str
        </path/to/alignment.sam>
    bam : str
        </path/to/output/alignment.bam>
    cpus : int, optional
        Number of processors to be used. By default uses all the processors of the system

    Raises
    ------
    TypeError
        cpus must be an integer greater than zero
    FileNotFoundError
        Specified path is incorrect or the file is empty
    subprocess.CalledProcessError
        Samtools did not run successfully
    """

    if type(cpus) is not int or cpus <= 0:
        raise TypeError(f"cpus must be an integer greater than zero! Given: {cpus}")
    if not os.path.exists(sam) or not os.path.getsize(sam):
        raise FileNotFoundError(
            f"The specified path: {sam} is either incorrect or the file is empty"
        )
    with tempfile.TemporaryDirectory() as tempdir:  # this will delete the temporary files even if program stops in between
        view_fname = "view_" + os.path.basename(bam)
        view_fpath = os.path.join(tempdir, view_fname)
        sort_fpath = os.path.join(tempdir, os.path.basename(bam))
        # See https://stackoverflow.com/questions/13332268/how-to-use-subprocess-command-with-pipes
        cmd_view = [
            "samtools",
            "view",
            "-@",
            str(cpus),
            "-bS",
            sam,
            "-o",
            view_fpath,
        ]
        cmd_sort = [
            "samtools",
            "sort",
            "-@",
            str(cpus),
            view_fpath,
            "-o",
            sort_fpath,
        ]
        logger.debug(" ".join(cmd_view))
        logger.debug(" ".join(cmd_sort))
        subprocess.run(
            cmd_view, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
        )
        subprocess.run(
            cmd_sort, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True,
        )
        shutil.move(sort_fpath, bam)
    return bam


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="%(asctime)s : %(name)s : %(levelname)s : %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    parser = argparse.ArgumentParser(
        description="Takes a sam file, sorts it and returns the output to a bam file"
    )
    parser.add_argument("sam", help="</path/to/alignment.sam>", type=str)
    parser.add_argument("bam", help="</path/to/output/alignment.bam>", type=str)
    parser.add_argument(
        "--cpus", help="Number of processors to use", default=mp.cpu_count(), type=int
    )
    args = parser.parse_args()
    sort(args.sam, args.bam, args.cpus)


if __name__ == "__main__":

    main()
