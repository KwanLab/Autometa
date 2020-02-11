#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script containing wrapper functions for bedtools.
"""


import logging
import os
import subprocess

import pandas as pd

logger = logging.getLogger(__name__)


def genomecov(ibam, lengths, out):
    """Run bedtools genomecov with input `ibam` and `lengths` to retrieve
    metagenome coverages.

    Parameters
    ----------
    ibam : str
        </path/to/indexed/BAM/file.ibam>. Note: BAM *must* be sorted by position.
    lengths : str
        </path/to/genome/lengths.tsv> tab-delimited cols=[contig,length]
    out : str
        </path/to/out.bed>
        The bedtools genomecov output is a tab-delimited file with the following columns:
        1. Chromosome
        2. Depth of coverage
        3. Number of bases on chromosome with that coverage
        4. Size of chromosome
        5. Fraction of bases on that chromosome with that coverage
        See also: http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    FileExistsError
        `out` file already exists
    OSError
        Why the exception is raised.

    """
    cmd = f'bedtools genomecov -ibam {ibam} -g {lengths}'
    if os.path.exists(out):
        raise FileExistsError(out)
    with open(os.devnull,'w') as stderr, open(out,'w') as stdout:
        retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
    if retcode or not os.path.exists(out) or os.stat(out).st_size == 0:
        raise OSError(f'bedtools failed: {cmd}')

def parse(bed, out=None):
    """Calculate coverages from bed file.

    Parameters
    ----------
    bed : str
        </path/to/file.bed>
    out : str
        if provided will write to `out`. I.e. </path/to/coverage.tsv>

    Returns
    -------
    pd.DataFrame
        index='contig', col='coverage'

    Raises
    -------
    FileNotFoundError
        `bed` does not exist

    """
    if out and os.path.exists(out):
        cols = ['contig','coverage']
        return pd.read_csv(out, sep='\t', usecols=cols, index_col='contig')
    if not os.path.exists(bed):
        raise FileNotFoundError(bed)
    names = ['contig','depth','bases','length','breadth']
    df = pd.read_csv(bed, sep='\t', names=names, index_col='contig')
    criterion1 = df.depth != 0
    criterion2 = df.index != 'genome'
    df = df[criterion1 & criterion2]
    df = df.assign(total_breadth=lambda x: x.depth * x.bases)
    dff = df.groupby('contig')['total_breadth', 'bases'].sum()
    dff = dff.assign(coverage=lambda x: x.total_breadth/x.bases)
    if out:
        dff.to_csv(out, sep='\t', index=True, header=True)
        logger.debug(f'{out} written')
    logger.debug(f'{os.path.basename(out)} shape: {dff.shape}')
    return dff['coverage']

def main(args):
    genomecov(ibam=args.ibam, lengths=args.lengths, out=args.bed)
    df = parse(bed=args.bed, out=args.coverage)

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    parser = argparse.ArgumentParser()
    parser.add_argument('ibam', help='</path/to/BAM/alignment.bam>')
    parser.add_argument(
        'lengths',
        help='</path/to/genome/lengths.tsv> tab-delimited cols=[contig,length]')
    parser.add_argument('bed',
        help='</path/to/alignment.bed> tab-delimited cols=[contig,length]')
    parser.add_argument('--coverage', help='</path/to/coverage.tsv>')
    args = parser.parse_args()
    main(args)
