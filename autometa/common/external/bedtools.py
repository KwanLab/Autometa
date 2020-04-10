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
Script containing wrapper functions for bedtools.
"""


import logging
import os
import subprocess

import pandas as pd

logger = logging.getLogger(__name__)


def genomecov(ibam, lengths, out, force=False):
    """Run bedtools genomecov with input `ibam` and `lengths` to retrieve
    metagenome coverages.

    Parameters
    ----------
    ibam : str
        </path/to/indexed/BAM/file.ibam>. Note: BAM *must* be sorted by position.
    lengths : str
        </path/to/genome/lengths.tsv> tab-delimited cols=[contig,length]
    out : str
        </path/to/alignment.bed>
        The bedtools genomecov output is a tab-delimited file with the following columns:
        1. Chromosome
        2. Depth of coverage
        3. Number of bases on chromosome with that coverage
        4. Size of chromosome
        5. Fraction of bases on that chromosome with that coverage
        See also: http://bedtools.readthedocs.org/en/latest/content/tools/genomecov.html
    force : bool
        force overwrite of `out` if it already exists (default is False).

    Returns
    -------
    str
        </path/to/alignment.bed>

    Raises
    -------
    FileExistsError
        `out` file already exists and force is False
    OSError
        Why the exception is raised.

    """
    cmd = f'bedtools genomecov -ibam {ibam} -g {lengths}'
    if os.path.exists(out) and not force:
        logger.debug(f'{out} already exists. skipping...')
        return out
    with open(os.devnull,'w') as stderr, open(out,'w') as stdout:
        retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
    if retcode or not os.path.exists(out) or os.stat(out).st_size == 0:
        raise ChildProcessError(f'bedtools failed: {cmd}')
    return out

def parse(bed, out=None, force=False):
    """Calculate coverages from bed file.

    Parameters
    ----------
    bed : str
        </path/to/file.bed>
    out : str
        if provided will write to `out`. I.e. </path/to/coverage.tsv>
    force : bool
        force overwrite of `out` if it already exists (default is False).

    Returns
    -------
    pd.DataFrame
        index='contig', col='coverage'

    Raises
    -------
    ValueError
        `out` incorrectly formatted to be read as pandas DataFrame.
    FileNotFoundError
        `bed` does not exist

    """
    if out and os.path.exists(out) and not os.stat(out).st_size == 0:
        try:
            cols = ['contig','coverage']
            return pd.read_csv(out, sep='\t', usecols=cols, index_col='contig')
        except ValueError as err:
            raise ValueError(f'InvalidTableFormat: {out}')
    if not os.path.exists(bed):
        raise FileNotFoundError(bed)
    names = ['contig','depth','bases','length','depth_fraction']
    df = pd.read_csv(bed, sep='\t', names=names, index_col='contig')
    criterion1 = df.depth != 0
    criterion2 = df.index != 'genome'
    df = df[criterion1 & criterion2]
    df = df.assign(depth_product=lambda x: x.depth * x.bases)
    dff = df.groupby('contig')['depth_product', 'bases'].sum()
    dff = dff.assign(coverage=lambda x: x.depth_product/x.bases)
    if out and (not os.path.exists(out) or (os.path.exists(out) and force)):
        dff.to_csv(out, sep='\t', index=True, header=True)
        logger.debug(f'{out} written')
    logger.debug(f'{os.path.basename(out)} shape: {dff.shape}')
    return dff[['coverage']]

def main(args):
    bed = genomecov(ibam=args.ibam, lengths=args.lengths, out=args.bed, force=args.force_bed)
    df = parse(bed=bed, out=args.coverage, force=args.force_cov)

if __name__ == '__main__':
    #start_parsing
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
    parser.add_argument('--force-bed', help='force overwrite `bed`',
        action='store_true',default=False)
    parser.add_argument('--force-cov', help='force overwrite `--coverage`',
        action='store_true',default=False)
    args = parser.parse_args()
    #end_parsing
    main(args)
