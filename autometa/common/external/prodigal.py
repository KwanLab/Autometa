#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to retrieve orfs from provided assembly using prodigal
"""

import gzip
import logging
import os
import subprocess
import shutil

from glob import glob


logger = logging.getLogger(__name__)


def run(assembly, nucls_out, prots_out, force=False,cpus=0,parallel=True):
    """Calls ORFs from provided input assembly

    Parameters
    ----------
    assembly : str
        </path/to/assembly.fasta>
    nucls_out : str
        </path/to/nucls.out>
    prots_out : str
        </path/to/prots.out>
    force : bool
        overwrite outfpath if it already exists (the default is False).
    cpus : int
        num `cpus` to use. By default will run as many `cpus` as possible
    parallel : bool
        Will parallelize prodigal using GNU parallel (the default is True).

    Returns
    -------
    2-Tuple
        (`nucls_out`, `prots_out`)

    Raises
    -------
    OSError
        Prodigal Failed
    """
    if not os.path.exists(assembly):
        raise FileNotFoundError(f'{assembly} does not exists!')
    if assembly.endswith('.gz'):
        with gzip.open(assembly) as fh:
            lines = ''
            for line in fh:
                lines += line.decode()
        assembly = assembly.rstrip('.gz')
        with open(assembly,'w') as fh:
            fh.write(lines)
    for fpath in [nucls_out, prots_out]:
        if os.path.exists(nucls_out) and not force:
            raise FileExistsError(f'{fpath} To overwrite use --force')
    if parallel:
        outdir = os.path.dirname(os.path.realpath(nucls_out))
        # parallel log should indicate time & dataset. i.e. time_dataset.prodigal.parallel.log
        log = os.path.join(outdir, 'prodigal.parallel.log')
        outprefix = os.path.splitext(os.path.basename(nucls_out))[0]
        tmpdir = os.path.join(outdir,'tmp')
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
        tmpnucl = '.'.join([outprefix, '{#}', 'fna'])
        tmpprot = '.'.join([outprefix, '{#}', 'faa'])
        tmpnucl_fpath = os.path.join(tmpdir, tmpnucl)
        tmpprot_fpath = os.path.join(tmpdir, tmpprot)
        jobs = f'-j{cpus}'
        cmd = [
            'parallel',
            '--retries',
            '4',
            '--joblog',
            log,
            jobs,
            '--pipe',
            '--recstart',
            '\'>\'',
            '--linebuffer',
            'prodigal',
            '-a',tmpprot_fpath,
            '-d',tmpnucl_fpath,
            '-q',
            '-p','meta',
            '-o',os.devnull,
            '<',assembly,
        ]
    else:
        cmd = [
            'prodigal',
            '-i',assembly,
            '-a',prots_out,
            '-d',nucls_out,
            '-p','meta',
            '-m',
            '-q',
        ]
    cmd = [str(arg) for arg in cmd]
    logger.debug(f'cmd: {" ".join(cmd)}')
    if parallel:
        returncode = subprocess.call(" ".join(cmd), shell=True)
        tmpfpaths = glob(os.path.join(tmpdir,'*.faa'))
        lines = ''
        for fp in tmpfpaths:
            with open(fp) as fh:
                for line in fh:
                    lines += line
        out = open(prots_out, 'w')
        out.write(lines)
        out.close()
        tmpfpaths = glob(os.path.join(tmpdir, '*.fna'))
        lines = ''
        for fp in tmpfpaths:
            with open(fp) as fh:
                for line in fh:
                    lines += line
        out = open(nucls_out, 'w')
        out.write(lines)
        out.close()
        shutil.rmtree(tmpdir)
    else:
        with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
            proc = subprocess.run(cmd, stdout=stdout, stderr=stderr)
            returncode = proc.returncode
    if returncode:
        logger.warning(f'Args:{cmd} ReturnCode:{returncode}')
    for fp in [nucls_out, prots_out]:
        if not os.path.exists(fp):
            raise OSError(f'{fp} not written')
    return nucls_out, prots_out

def main(args):
    if args.verbose:
        logger.setLevel(logger.DEBUG)
    nucls_out, prots_out = run(
        assembly=args.assembly,
        nucls_out=args.nucls_out,
        prots_out=args.prots_out,
        force=args.force,
        cpus=args.cpus,
        parallel=args.noparallel)
    logger.info(f'written:\nnucls fpath: {nucls_out}\nprots fpath: {prots_out}')

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Calls ORFs with provided input assembly')
    parser.add_argument('assembly', help='</path/to/assembly>')
    parser.add_argument('nucls_out', help='</path/to/nucls.out>')
    parser.add_argument('prots_out', help='</path/to/prots.out>')
    parser.add_argument('--force', help="force overwrite of ORFs out filepaths",
        action='store_true')
    parser.add_argument('--cpus', help='num cpus to use', default=0)
    parser.add_argument('--noparallel', help="Disable GNU parallel",
        action='store_false', default=True)
    parser.add_argument('--verbose', help="add verbosity", action='store_true')
    args = parser.parse_args()
    main(args)
