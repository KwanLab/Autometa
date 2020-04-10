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

Functions to retrieve orfs from provided assembly using prodigal
"""

import gzip
import logging
import os
import subprocess
import shutil

from glob import glob
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from autometa.config.environ import get_versions


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
        num `cpus` to use. **Default (cpus=0) will run as many `cpus` as possible**
    parallel : bool
        Will parallelize prodigal using GNU parallel (the default is True).

    Returns
    -------
    2-Tuple
        (`nucls_out`, `prots_out`)

    Raises
    -------
    FileExistsError
        `nucls_out` or `prots_out` already exists
    ChildProcessError
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
        if os.path.exists(fpath) and not force:
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
            '--retries','4',
            '--joblog',log,
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
            '-m',
            '-o',os.devnull,
            '<',assembly,
            '2>',os.devnull,
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
        try:
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
        except Exception as err:
            # COMBAK: Should probably be more descriptive as to what errors could occur here.
            logger.exception(err)
        finally:
            shutil.rmtree(tmpdir)
    else:
        with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
            proc = subprocess.run(cmd, stdout=stdout, stderr=stderr)
            returncode = proc.returncode
    if returncode:
        logger.warning(f'Args:{cmd} ReturnCode:{returncode}')
        # COMBAK: Check all possible return codes for GNU parallel
    for fp in [nucls_out, prots_out]:
        if not os.path.exists(fp) or os.stat(fp).st_size == 0:
            raise ChildProcessError(f'{fp} not written')
        try:
            with open(fp) as fh:
                for _ in SimpleFastaParser(fh):
                    pass
        except (IOError, ValueError):
            raise IOError(f'InvalidFileFormat: {fp}')
    return nucls_out, prots_out

def contigs_from_headers(fpath):
    """Get ORF id to contig id translations using prodigal assigned ID from
    description.

    First determines if all of ID=3495691_2 from description is in header.
    "3495691_2" represents the 3,495,691st gene in the 2nd sequence.

    Example
    -------
    .. code-block:: python

        #: prodigal versions < 2.6 record
        >>>record.id
        'k119_1383959_3495691_2'

        >>>record.description
        'k119_1383959_3495691_2 # 688 # 1446 # 1 # ID=3495691_2;partial=01;start_type=ATG;rbs_motif=None;rbs_spacer=None'

        >>>record.description.split('#')[-1].split(';')[0].strip()
        'ID=3495691_2'

        >>>orf_id = '3495691_2'
        '3495691_2'

        >>>record.id.replace(f'_{orf_id}', '')
        'k119_1383959'

        #: prodigal versions >= 2.6 record
        >>>record.id
        'k119_1383959_2'
        >>>record.id.rsplit('_',1)[0]
        'k119_1383959'

    Parameters
    ----------
    fpath : str
        </path/to/prodigal/called/orfs.fasta>

    Returns
    -------
    dict
        contigs translated from prodigal ORF description.  {orf_id:contig_id, ...}

    """
    version = get_versions('prodigal')
    if version.count('.') >= 2:
        version = float('.'.join(version.split('.')[:2]))
    else:
        version = float(version)
    translations = {}
    for record in SeqIO.parse(fpath, 'fasta'):
        if version < 2.6:
            orf_id = record.description.split('#')[-1].split(';')[0].strip().replace('ID=','')
            contig_id = record.id.replace(f'_{orf_id}', '')
        else:
            contig_id = record.id.rsplit('_',1)[0]
        translations.update({record.id:contig_id})
    return translations

def orf_records_from_contigs(contigs, fpath):
    """Retrieve list of *ORFs headers* from `contigs`. Prodigal annotated ORFs
    are required as the input `fpath`.

    Parameters
    ----------
    contigs: iterable
        iterable of contigs from which to retrieve ORFs
    fpath : str
        </path/to/prodigal/called/orfs.fasta>

    Returns
    -------
    list
        ORF SeqIO.SeqRecords from provided `contigs`. i.e. [SeqRecord, ...]

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    version = get_versions('prodigal')
    if version.count('.') >= 2:
        version = float('.'.join(version.split('.')[:2]))
    else:
        version = float(version)

    records = []
    for record in SeqIO.parse(fpath, 'fasta'):
        if version < 2.6:
            orf_id = record.description.split('#')[-1].split(';')[0].strip().replace('ID=','')
            contig_id = record.id.replace(f'_{orf_id}', '')
        else:
            contig_id = record.id.rsplit('_',1)[0]
        if contig_id not in contigs:
            continue
        records.append(record)
    return records

def main(args):
    if args.verbose:
        logger.setLevel(logger.DEBUG)
    nucls_out, prots_out = run(
        assembly=args.assembly,
        nucls_out=args.nucls_out,
        prots_out=args.prots_out,
        force=args.force,
        cpus=args.cpus,
        parallel=args.parallel)
    logger.info(f'written:\nnucls fpath: {nucls_out}\nprots fpath: {prots_out}')

if __name__ == '__main__':
    #start_parsing
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser(usage = 'prodigal.py' ,
        description = 'Calls ORFs with provided input assembly')
    parser.add_argument('assembly', help='</path/to/assembly>', type=str)
    parser.add_argument('nucls_out', help='</path/to/nucls.out>', type=str)
    parser.add_argument('prots_out', help='</path/to/prots.out>', type=str)
    parser.add_argument('--force', help="force overwrite of ORFs out filepaths",
        action='store_true')
    parser.add_argument('--cpus', help='num cpus to use', type=int, default=0)
    parser.add_argument('--parallel', help="Enable GNU parallel",
        action='store_true', default=False)
    parser.add_argument('--verbose', help="add verbosity", action='store_true')
    args = parser.parse_args()
    #end_parsing
    main(args)
