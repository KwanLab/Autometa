#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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

Script containing wrapper functions for samtools
"""


import logging
import os
import subprocess

import pandas as pd
import multiprocessing as mp

logger = logging.getLogger(__name__)


def run(cmd):
    """Run `cmd` via subprocess.

    Parameters
    ----------
    cmd : str
        Executable input str

    Returns
    -------
    bool
        True if retcode == 0. i.e. Success.
        False if retcode != 0. i.e. Failed.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    logger.debug(f'cmd: {cmd}')
    with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
        retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
    if not retcode:
        return True
    else:
        logger.warning(f'args:{cmd} retcode:{retcode}')
        return False

def sort(sam, out, nproc=mp.cpu_count()):
    """
    Views the sam file and then sorts the alighnments by leftmost coordinates.
    
    Parameters
    ----------
    sam : str
        </path/to/alignment.sam>
    out : str
        </path/to/output/file
    nproc : int, optional
        Number of processors to be used. By default uses all the processors of the system. 

    """
    cmd = f'samtools view -@{nproc} -bS {sam} | samtools sort -@{nproc} -o {out}'
    run(cmd)

def main(args):
    sort(args.sam, args.out, args.nproc)

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    parser = argparse.ArgumentParser()
    parser.add_argument('--sam', help='</path/to/alignment.sam>')
    parser.add_argument('--out', help='</path/to/output/file')
    parser.add_argument('--nproc', help='Number of processors to use', default=mp.cpu_count(), type=int)
    args = parser.parse_args()
    main(args)
