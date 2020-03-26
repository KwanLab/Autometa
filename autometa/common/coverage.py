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
Autometa Coverage
"""


import logging
import os
import shutil
import tempfile

import pandas as pd

from Bio import SeqIO

from autometa.common.external import bowtie
from autometa.common.external import samtools
from autometa.common.external import bedtools

logger = logging.getLogger(__name__)


def from_spades_names(records):
    """Retrieve coverages from SPAdes scaffolds headers.

    Example SPAdes header : NODE_83_length_162517_cov_224.639

    Parameters
    ----------
    records : list
        [SeqRecord,...]

    Returns
    -------
    pd.Series
        index=contig, name='coverage', dtype=float

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    return pd.Series(
        {record.id:record.id.split('_cov_')[-1] for record in records},
        name='coverage',
        dtype=float)

def make_length_table(fasta, out):
    seqs = {r.id:len(r) for r in SeqIO.parse(fasta, 'fasta')}
    length_s = pd.Series(seqs, name='length')
    length_s.index.name = 'contig'
    length_s.to_csv(out, sep='\t', index=True, header=True)
    return out

def get(fasta, out, fwd_reads=None, rev_reads=None, sam=None, bam=None, lengths=None,
    bed=None, nproc=1):
    """Get coverages for assembly `fasta` file using provided files:

    Either:
        `fwd_reads` and `rev_reads`
    or:
        `sam`
    or:
        `bam`
    or:
        `bed`

    Will begin coverage calculation based on files provided checking in the
    following order:
        1. `bed`
        2. `bam`
        3. `sam`
        4. `fwd_reads` and `rev_reads`


    Parameters
    ----------
    fasta : str
        </path/to/assembly.fasta>
    out : str
        </path/to/output/coverage.tsv>
    fwd_reads : str
        </path/to/fwd_reads.fastq>
    rev_reads : str
        </path/to/rev_reads.fastq>
    sam : str
        </path/to/alignments.sam>
    bam : str
        </path/to/alignments.bam>
    lengths : str
        </path/to/lengths.tsv>
    bed : str
        </path/to/alignments.bed>

    Returns
    -------
    pd.DataFrame
        index=contig cols=['coverage']
    """
    if os.path.exists(out) and os.stat(out).st_size > 0:
        # COMBAK: checksum comparison [checkpoint]
        logger.debug(f'coverage ({out}) already exists. skipping...')
        cols = ['contig','coverage']
        return pd.read_csv(out, sep='\t', usecols=cols, index_col='contig')
    try:
        outdir = os.path.dirname(out)
        tempdir = tempfile.mkdtemp(suffix=None, prefix='cov-alignments', dir=outdir)
        bed = bed if bed else os.path.join(tempdir, 'alignment.bed')
        bam = bam if bam else os.path.join(tempdir, 'alignment.bam')
        lengths = lengths if lengths else os.path.join(tempdir, 'lengths.tsv')
        sam = sam if sam else os.path.join(tempdir, 'alignment.sam')
        db = os.path.join(tempdir, 'alignment.db')
        if os.path.exists(bed):
            return bedtools.parse(bed, out)
        if os.path.exists(bam):
            if not os.path.exists(lengths):
                lengths = make_length_table(fasta, lengths)
            bedtools.genomecov(bam, lengths, bed)
            return bedtools.parse(bed, out)
        if os.path.exists(sam):
            samtools.sort(sam, bam, nproc=nproc)
            if not os.path.exists(lengths):
                lengths = make_length_table(fasta, lengths)
            bedtools.genomecov(bam, lengths, bed)
            return bedtools.parse(bed, out)
        if not fwd_reads or not rev_reads:
            raise ValueError(f'{fwd_reads} and {rev_reads} are required if no other alignments are specified!')
        bowtie.build(fasta, db)
        bowtie.align(db, sam, fwd_reads, rev_reads, nproc=nproc)
        samtools.sort(sam, bam, nproc=nproc)
        if not os.path.exists(lengths):
            lengths = make_length_table(fasta, lengths)
        bedtools.genomecov(bam, lengths, bed)
        return bedtools.parse(bed, out)
    except Exception as err:
        logger.exception(err)
    finally:
        shutil.rmtree(tempdir, ignore_errors=True)

def main(args):
    get(fasta=args.assembly,
        fwd_reads=args.fwd_reads,
        rev_reads=args.rev_reads,
        sam=args.sam,
        bam=args.bam,
        lengths=args.lengths,
        bed=args.bed,
        nproc=args.nproc,
        out=args.out)

if __name__ == '__main__':
    import argparse
    import logging as logger
    import multiprocessing as mp
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Autometa Coverage')
    parser.add_argument('-f','--assembly', help='</path/to/metagenome.fasta>', required=True)
    parser.add_argument('-1', '--fwd-reads', help='</path/to/forwards-reads.fastq>')
    parser.add_argument('-2', '--rev-reads', help='</path/to/reverse-reads.fastq>')
    parser.add_argument('--sam', help='</path/to/alignments.sam>')
    parser.add_argument('--bam', help='</path/to/alignments.bam>')
    parser.add_argument('--lengths', help='</path/to/lengths.tsv>')
    parser.add_argument('--bed', help='</path/to/alignments.bed>')
    parser.add_argument('--nproc', help=f'Num processors to use. (default: {mp.cpu_count()})', default=mp.cpu_count())
    parser.add_argument('--out', help='</path/to/coverages.tsv>')
    args = parser.parse_args()
    main(args)
