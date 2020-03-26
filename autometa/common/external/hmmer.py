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
Functions related to running hmmer on metagenome sequences
"""


import logging
import os
import subprocess
import shutil

import pandas as pd

from glob import glob
from Bio import SeqIO

from autometa.common.external import prodigal


logger = logging.getLogger(__name__)


def hmmscan(orfs, hmmdb, outfpath, cpus=0, force=False, parallel=True, log=None):
    """Runs hmmscan on dataset ORFs and provided hmm database.

    Parameters
    ----------
    orfs : str
        </path/to/orfs.faa>
    hmmdb : str
        </path/to/hmmpressed/database.hmm>
    outfpath : str
        </path/to/output.hmmscan.tsv>
    cpus : int
        Num. cpus to use. 0 will run as many cpus as possible (the default is 0).
    force : bool
        Overwrite existing `outfpath` (the default is False).
    parallel : bool
        Will parallelize hmmscan using GNU parallel (the default is True).
    log : str
        </path/to/parallel.log> (the default is None). If provided will write
        parallel log to `log`.

    Returns
    -------
    str
        </path/to/output.hmmscan.tsv>

    Raises
    -------
    FileExistsError
        `outfpath` already exists
    OSError
        hmmscan failed
    """
    # OPTIMIZE: we want to extend parallel to grid computing (workqueue?) via --sshlogin?
    if os.path.exists(outfpath) and os.stat(outfpath).st_size > 0 and not force:
        raise FileExistsError(f'{outfpath}. Use force to overwrite!')
    if parallel:
        outdir = os.path.dirname(os.path.realpath(outfpath))
        outprefix = os.path.splitext(os.path.basename(outfpath))[0]
        tmpdir = os.path.join(outdir,'tmp')
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
        tmpfname = '.'.join([outprefix, '{#}', 'txt'])
        tmpfpath = os.path.join(tmpdir, tmpfname)
        jobs = f'-j{cpus}'
        cmd = [
            'parallel',
            '--retries',
            '4',
            jobs,
            '--linebuffer',
            '--pipe',
            '--recstart',
            '\'>\'',
            'hmmscan',
            '-o',os.devnull,
            '--tblout',
            tmpfpath,
            hmmdb,
            '-',
            '<',orfs,
            '2>',os.devnull,
        ]
        if log:
            cmd.insert(3,log)
            cmd.insert(3,'--joblog')
    else:
        cmd = ['hmmscan', '--cpu', str(cpus), '--tblout', outfpath, hmmdb, orfs]
    logger.debug(f'cmd: {" ".join(cmd)}')
    if parallel:
        returncode = subprocess.call(' '.join(cmd), shell=True)
        tmpfpaths = glob(os.path.join(tmpdir, '*.txt'))
        lines = ''
        for fp in tmpfpaths:
            with open(fp) as fh:
                for line in fh:
                    lines += line
        out = open(outfpath, 'w')
        out.write(lines)
        out.close()
        shutil.rmtree(tmpdir)
    else:
        with open(os.devnull, 'w') as STDOUT, open(os.devnull, 'w') as STDERR:
            proc = subprocess.run(cmd, stdout=STDOUT, stderr=STDERR)
        returncode = proc.returncode
    if returncode == 141:
        logger.warning(f'Make sure your hmm profiles are pressed!')
        logger.warning(f'hmmpress -f {hmmdb}')
        logger.error(f'Args:{cmd} ReturnCode:{returncode}')
        raise OSError(f'hmmscan: Args:{cmd} ReturnCode:{returncode}')
    if not os.path.exists(outfpath):
        raise OSError(f'{outfpath} not written.')
    return outfpath

def filter_markers(infpath, outfpath, cutoffs, prodigal_annotations=None, force=False):
    """Filter markers from hmmscan output table that are above cutoff values.

    Parameters
    ----------
    infpath : str
        </path/to/hmmscan.tsv>
    outfpath : str
        </path/to/output.markers.tsv>
    cutoffs : str
        </path/to/cutoffs.tsv>
    prodigal_annotations : str
        Default will attempt to translate recovered qseqids to contigs
        </path/to/prodigal/called/orfs.fasta>
    force : bool
        Overwrite existing `outfpath` (the default is False).

    Returns
    -------
    str
        </path/to/output.markers.tsv>

    Raises
    -------
    FileNotFoundError
        `infpath` or `cutoffs` not found
    FileExistsError
        `outfpath` already exists and force=False
    AssertionError
        No returned markers pass the cutoff thresholds. I.e. final df is empty.
    """
    for fp in [infpath, cutoffs]:
        if not os.path.exists(fp):
            raise FileNotFoundError(f'{fp} not found')
    if os.path.exists(outfpath) and os.stat(outfpath).st_size > 0 and not force:
        raise FileExistsError(f'{outfpath} already exists')
    hmmtab_header = ['sname','sacc','orf','score']
    col_indices = [0,1,2,5]
    columns = {i:k for i,k in zip(col_indices,hmmtab_header)}
    df = pd.read_csv(
        infpath,
        sep='\s+',
        usecols=col_indices,
        header=None,
        comment='#',
    )
    df.rename(columns=columns, inplace=True)
    # NaN may result from parsing issues while merging parallel results
    df.dropna(inplace=True)
    df['cleaned_sacc'] = df['sacc'].map(lambda acc: acc.split('.')[0])
    dff = pd.read_csv(cutoffs, sep='\t', index_col='accession')
    mdf = pd.merge(df,dff,how='left',left_on='cleaned_sacc',right_on='accession')
    mdf = mdf[mdf['score'] >= mdf['cutoff']]
    if mdf.empty:
        raise AssertionError(f'No markers in {infpath} pass cutoff thresholds')
    cols = ['orf','sacc','sname','score','cutoff']
    mdf = mdf[cols]
    if prodigal_annotations:
        logger.debug('Retrieving ORF->contig translations from ORF Caller')
        translations = prodigal.contigs_from_headers(prodigal_annotations)
        translater = lambda x: translations.get(x, x.rsplit('_',1)[0])
    else:
        translater = lambda x: x.rsplit('_',1)[0]
    mdf['contig'] = mdf['orf'].map(translater)
    mdf.set_index('contig', inplace=True)
    mdf.to_csv(outfpath, sep='\t', index=True, header=True)
    return outfpath

def main(args):
    if args.verbose:
        logger.setLevel(logger.DEBUG)

    try:
        result = hmmscan(
            orfs=args.orfs,
            hmmdb=args.hmmdb,
            outfpath=args.hmmscan,
            cpus=args.cpus,
            force=args.force,
            parallel=args.noparallel,
            log=args.log)
    except FileExistsError as err:
        logger.debug(err)
        result = args.hmmscan

    result = filter_markers(
        infpath=result,
        outfpath=args.markers,
        cutoffs=args.cutoffs,
        prodigal_annotations=args.orfs,
        force=args.force)


if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    parser = argparse.ArgumentParser('Retrieves markers with provided input assembly')
    parser.add_argument('orfs', help='</path/to/assembly.orfs.faa>')
    parser.add_argument('hmmdb', help='</path/to/hmmpressed/hmmdb>')
    parser.add_argument('cutoffs', help='</path/to/hmm/cutoffs.tsv>')
    parser.add_argument('hmmscan', help='</path/to/hmmscan.out>')
    parser.add_argument('markers', help='</path/to/markers.tsv>')
    parser.add_argument('--log', help='</path/to/parallel.log>')
    parser.add_argument('--force', help="force overwrite of out filepath",
        action='store_true')
    parser.add_argument('--cpus', help='num cpus to use',default=0)
    parser.add_argument('--noparallel',help="Disable GNU parallel", action='store_false')
    parser.add_argument('--verbose', help="add verbosity", action='store_true')
    args = parser.parse_args()
    main(args)
