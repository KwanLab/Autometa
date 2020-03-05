#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script containing wrapper functions for samtools
"""


import logging
import os
import subprocess

import pandas as pd

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
    if retcode:
        logger.warning(f'args:{cmd} retcode:{retcode}')
        return False
    return True

def sort(sam, out, nproc=0):
    cmd = f'samtools view -@{nproc} -bS {sam} | samtools sort -@{nproc} -o {out}'
    run(cmd)

def depth(bam, records):
    cmds = [f'samtools depth -r {record.id} {bam}' for record in records]
    raise NotImplementedError

def main(args):
    sort(args.sam)
    depth(args.bam, args.seqs)

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam')
    parser.add_argument('--sam')
    parser.add_argument('--seqs')
    args = parser.parse_args()
    main(args)
