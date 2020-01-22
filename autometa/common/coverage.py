#!/usr/bin/env python
"""
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

def get(fasta, fwd_reads, rev_reads, nproc=1, out=None):
    """Get coverages for assembly `fasta` file using provided `fwd_reads` and
    `rev_reads`.

    Parameters
    ----------
    fasta : str
        </path/to/assembly.fasta>
    fwd_reads : str
        </path/to/fwd_reads.fastq>
    rev_reads : str
        </path/to/rev_reads.fastq>
    out : str
        </path/to/output/coverage.tsv>

    Returns
    -------
    pd.DataFrame
        index=contig cols=['coverage']
    """
    try:
        outdir = os.path.dirname(out)
        tempdir = tempfile.mkdtemp(suffix=None, prefix='cov-alignments', dir=outdir)
        db = os.path.join(tempdir, 'alignment.db')
        sam = os.path.join(tempdir, 'alignment.sam')
        bam = os.path.join(tempdir, 'alignment.bam')
        bed = os.path.join(tempdir, 'alignment.bed')
        lengths = os.path.join(tempdir, 'lengths.tsv')
        bowtie.build(fasta, db)
        bowtie.align(db, sam, fwd_reads, rev_reads, nproc=nproc)
        samtools.sort(sam, bam, nproc=nproc)
        # samtools.depth(bam, records)
        seqs = {r.id:len(r) for r in SeqIO.parse(fasta, 'fasta')}
        length_s = pd.Series(seqs, name='length')
        length_s.index.name = 'contig'
        length_s.to_csv(lengths, sep='\t', index=True, header=True)
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
    parser.add_argument('-1', '--fwd-reads', help='</path/to/forwards-reads.fastq>',
        required=True)
    parser.add_argument('-2', '--rev-reads', help='</path/to/reverse-reads.fastq>',
        required=True)
    parser.add_argument('--nproc', help='Num processors to use.', default=mp.cpu_count())
    parser.add_argument('--out', help='</path/to/read/alignments.bam>')
    args = parser.parse_args()
    main(args)
