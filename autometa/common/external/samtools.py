#!/usr/bin/env python
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
    
    Raises
    ------
    TypeError
        nproc must be an integer greater than zero
    ChildProcessError
        Samtools did not run successfully, returns the return code from subprocessing call
    """

    if type(nproc) is not int or nproc <= 0:
        raise TypeError(f'nproc must be an integer greater than zero! Given: {nproc}')
    samtools_out_dir = os.path.dirname(os.path.abspath(bam))
    samtools_err = os.path.join(samtools_out_dir, 'samtools.err')
    samtools_out = os.path.join(samtools_out_dir, 'samtools.out')
    outfile = os.path.basename(bam)
    samtools_outfile = os.path.join(samtools_out_dir, outfile)
    with tempfile.TemporaryDirectory(suffix = None, prefix = 'samtools', dir = samtools_out_dir) as tempdir:    
        temp_outfile = os.path.join(tempdir, outfile)
        cmd = f'samtools view -@{nproc} -bS {sam} | samtools sort -@{nproc} -o {temp_outfile}'
        logger.debug(f'cmd: {cmd}')
        with open(samtools_out, 'w') as stdout, open(samtools_err, 'w') as stderr:
            retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
            if retcode != 0:
                raise ChildProcessError(f'Sort failed, retcode: {retcode}')
            else:
                shutil.move(temp_outfile, samtools_outfile)
        os.remove(samtools_err)
        os.remove(samtools_out)
          
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
