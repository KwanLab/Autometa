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
    FileNotFoundError
        Specified path is incorrect or the file is empty        
    ChildProcessError
        Samtools did not run successfully, returns the return code from subprocessing call
    """

    if type(nproc) is not int or nproc <= 0:
        raise TypeError(f'nproc must be an integer greater than zero! Given: {nproc}')
    if not os.path.exists(sam) or os.stat(sam).st_size == 0:
        raise FileNotFoundError(f'The specified path: {sam} is either incorrect or the file is empty')
    samtools_out_dir = os.path.dirname(os.path.abspath(bam))
    err = os.path.join(samtools_out_dir, 'samtools.err')
    out = os.path.join(samtools_out_dir, 'samtools.out')
    with tempfile.TemporaryDirectory() as tempdir:    
        temp_bam = os.path.join(tempdir, os.path.basename(bam))
        cmd = f'samtools view -@{nproc} -bS {sam} | samtools sort -@{nproc} -o {temp_bam}'
        logger.debug(f'cmd: {cmd}')
        with open(out, 'w') as stdout, open(err, 'w') as stderr:
            retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
        if retcode != 0:
            raise ChildProcessError(f'Sort failed, retcode: {retcode}')
        else:
            shutil.move(temp_bam, bam)
        os.remove(err)
        os.remove(out)
          
def main(args):
    sort(args.sam, args.bam, args.nproc)

if __name__ == '__main__':
    import argparse
    import logging as logger
    import multiprocessing as mp
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    parser = argparse.ArgumentParser(description = "Takes a sam file, sorts it and returns the output to a bam file")
    parser.add_argument('--sam', help='</path/to/alignment.sam>', type=str)
    parser.add_argument('--bam', help='</path/to/output/alignment.bam>', type=str)
    parser.add_argument('--nproc', help='Number of processors to use', default=mp.cpu_count(), type=int)
    args = parser.parse_args()
    main(args)
