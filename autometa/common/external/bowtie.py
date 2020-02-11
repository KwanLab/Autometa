#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Script containing wrapper functions for bowtie2.
"""


import logging
import os
import subprocess


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
        True if no returncode from subprocess.run else False

    """
    logger.debug(f'run: {cmd}')
    with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
        retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
    if retcode:
        logger.warning(f'Args:{cmd} ReturnCode:{retcode}')
        return False
    return True

def build(assembly, out):
    """Build bowtie2 index.

    Parameters
    ----------
    assembly : str
        </path/to/assembly.fasta>
    out : str
        </path/to/output/database>
        Note: Indices written will resemble </path/to/output/database.{#}.bt2>

    Returns
    -------
    str
        </path/to/output/database>

    Raises
    -------
    OSError
        bowtie2-build failed
    """
    cmd = f'bowtie2-build {assembly} {out}'
    success = run(cmd)
    if not success:
        raise OSError(f'{cmd} failed. {out} not written')
    return out

def align(db, sam, fwd_reads, rev_reads, nproc=0, **kwargs):
    """Align reads to bowtie2-index `db`.

    Parameters
    ----------
    db : str
        </path/to/prefix/bowtie2/database>. I.e. `db`.{#}.bt2
    sam : str
        </path/to/out.sam>
    fwd_reads : str
        </path/to/forward_reads.fastq>
    rev_reads : str
        </path/to/reverse_reads.fastq>
    nproc : int
        Num. processors to use (the default is 0).
    **kwargs : dict
        Additional optional args to supply to bowtie2. Must be in format:
        key = flag
        value = flag-value

    Returns
    -------
    str
        </path/to/out.sam>

    Raises
    -------
    OSError
        bowtie2 failed
    """
    exc = f'bowtie2 -x {db}'
    opts = '-q --phred33 --very-sensitive --no-unal'
    added_opts = [f'{k} {v}' for k,v in kwargs.items()] if kwargs else None
    nprocs = f'-p {nproc}'
    sam_out = f'-S {sam}'
    fwd_reads_param = f'-1 {fwd_reads}'
    rev_reads_param = f'-2 {rev_reads}'
    params = [exc,opts,nprocs,fwd_reads_param,rev_reads_param,sam_out]
    if added_opts:
        params += added_opts
    cmd =  ' '.join(params)
    success = run(cmd)
    if not success:
        raise OSError(f'{cmd} failed. {sam} not written')
    return sam

def main(args):
    db = build(args.assembly, args.database)
    sam = align(args.database, args.sam, args.fwd_reads, args.rev_reads, args.nproc)

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    parser = argparse.ArgumentParser()
    parser.add_argument('assembly', help='</path/to/assembly.fasta>')
    parser.add_argument('database', help='</path/to/alignment.database>')
    parser.add_argument('sam', help='</path/to/alignment.sam>')
    parser.add_argument('-1', '--fwd-reads', help='</path/to/forwards-reads.fastq>')
    parser.add_argument('-2', '--rev-reads', help='</path/to/reverse-reads.fastq>')
    parser.add_argument('--nproc', help='Num processors to use.', default=1)
    args = parser.parse_args()
    main(args)
