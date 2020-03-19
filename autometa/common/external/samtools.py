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
import shutil
import tempfile

import pandas as pd
import multiprocessing as mp

logger = logging.getLogger(__name__)


def run(cmd, bam):
    """
    Run `cmd` via subprocess.
    
    Parameters
    ----------
    cmd : str
        Executable input str
    bam : str
        </path/to/output/alignment.bam>
    
    Returns
    -------
    bool
        True if retcode == 0. i.e. function run successfully.
    
    Raises
    ------
    ChildProcessError
        Function did not run successfully, returns retcode.
    """
    
    logger.debug(f'cmd: {cmd}')
    log_samtools_dir = os.path.dirname(bam)
    tempdir = tempfile.mkdtemp(suffix=None, prefix='samtools', dir=log_samtools_dir)
    samtools_stderr=os.path.join(tempdir, 'samtools_stderr')
    samtools_stdout=os.path.join(tempdir, 'samtools_stdout')
    with open(samtools_stdout, 'w') as stdout, open(samtools_stderr, 'w') as stderr:
        retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
        if retcode != 0:
            raise ChildProcessError(f'Samtools not successfully run, retcode: {retcode}')
        shutil.rmtree(tempdir, ignore_errors=True)
        return True
    
def sort(sam, bam, nproc=mp.cpu_count()):
    """ 
    Views then sorts sam file by leftmost coordinates and outputs to bam.
    
    Parameters
    ----------
    sam : str
        </path/to/alignment.sam>
    bam : str
        </path/to/output/alignment.bam>
    nproc : int, optional
        Number of processors to be used. By default uses all the processors of the system 
    """
    
    cmd = f'samtools view -@{nproc} -bS {sam} | samtools sort -@{nproc} -o {bam}'
    run(cmd,bam)

def main(args):
    sort(args.sam, args.bam, args.nproc)

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    parser = argparse.ArgumentParser(description = "Takes a sam file, sorts it and returns the output to a bam file")
    parser.add_argument('--sam', help='</path/to/alignment.sam>', type=str)
    parser.add_argument('--bam', help='</path/to/output/alignment.bam>', type=str)
    parser.add_argument('--nproc', help='Number of processors to use', default=mp.cpu_count(), type=int)
    args = parser.parse_args()
    main(args)

