#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Calculates coverage of contigs
"""


import gzip
import logging
import os
import tempfile
from typing import List

import pandas as pd

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from autometa.common.external import bowtie
from autometa.common.external import samtools
from autometa.common.external import bedtools
from autometa.common import utilities


logger = logging.getLogger(__name__)


def from_spades_names(records: List[SeqRecord]) -> pd.DataFrame:
    """Retrieve coverages from SPAdes scaffolds headers.

    Example SPAdes header : NODE_83_length_162517_cov_224.639

    Parameters
    ----------
    records : List[SeqRecord]
        [SeqRecord,...]

    Returns
    -------
    pd.DataFrame
        index=contig, name='coverage', dtype=float

    """
    logger.info(f"Retrieving coverages from contig ID in {len(records):,} records")
    coverages = pd.Series(
        {record.id: record.id.split("_cov_")[-1] for record in records},
        name="coverage",
        dtype=float,
    )
    coverages.index.name = "contig"
    return coverages.to_frame()


@utilities.timeit
def get(
    fasta: str,
    out: str,
    from_spades: bool = False,
    fwd_reads: List[str] = None,
    rev_reads: List[str] = None,
    se_reads: List[str] = None,
    sam: str = None,
    bam: str = None,
    bed: str = None,
    cpus: int = 1,
) -> pd.DataFrame:
    """Get coverages for assembly `fasta` file using provided files or
    if the metagenome assembly was generated from SPAdes, use the k-mer coverages
    provided in each contig's header by specifying `from_spades=True`.

    Either `fwd_reads` and `rev_reads` and/or `se_reads` or,`sam`, or `bam`, or `bed`
    must be provided if `from_spades=False`.

    Note
    ----

        Will begin coverage calculation based on files provided checking in the
        following order:

        #. `bed`

        #. `bam`

        #. `sam`

        #. `fwd_reads` and `rev_reads` and `se_reads`

        Event sequence to calculate contig coverages:

        #. align reads to generate alignment.sam

        #. sort samfile to generate alignment.bam

        #. calculate assembly coverages to generate alignment.bed

        #. calculate contig coverages to generate coverage.tsv


    Parameters
    ----------

    fasta : str
        </path/to/assembly.fasta>

    out : str
        </path/to/output/coverages.tsv>

    from_spades : bool, optional
        If True, will attempt to parse record ids for coverage information.
        This is only compatible with SPAdes assemblies. (the Default is False).

    fwd_reads : List[str], optional
        [</path/to/forward_reads.fastq>, ...]

    rev_reads : List[str], optional
        [</path/to/reverse_reads.fastq>, ...]

    se_reads : List[str], optional
        [</path/to/single_end_reads.fastq>, ...]

    sam : str, optional
        </path/to/alignments.sam>

    bam : str, optional
        </path/to/alignments.bam>

    bed : str, optional
        </path/to/alignments.bed>

    cpus : int, optional
        Number of cpus to use for coverage calculation.

    Returns
    -------
    pd.DataFrame
        index=contig cols=['coverage']
    """

    if os.path.exists(out) and os.path.getsize(out):
        logger.debug(f"coverage ({out}) already exists. skipping...")
        cols = ["contig", "coverage"]
        return pd.read_csv(out, sep="\t", usecols=cols, index_col="contig")

    if from_spades:
        fh = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta)
        df = from_spades_names(records=[rec for rec in SeqIO.parse(fh, "fasta")])
        fh.close()
        df.to_csv(out, sep="\t", index=True, header=True)
        return df

    with tempfile.TemporaryDirectory() as tempdir:
        outdir = os.path.dirname(out)
        os.makedirs(outdir, exist_ok=True)
        tempdir = tempfile.mkdtemp(suffix=None, prefix="cov-alignments", dir=outdir)
        bed = bed if bed else os.path.join(tempdir, "alignment.bed")
        bam = bam if bam else os.path.join(tempdir, "alignment.bam")
        sam = sam if sam else os.path.join(tempdir, "alignment.sam")
        db = os.path.join(tempdir, "alignment.db")

        def parse_bed(bed=bed, out=out):
            return bedtools.parse(bed, out)

        def make_bed(bam=bam, bed=bed):
            bedtools.genomecov(bam, bed)

        def sort_samfile(sam=sam, bam=bam, cpus=cpus):
            samtools.sort(sam, bam, cpus=cpus)

        def align_reads(
            fasta=fasta,
            db=db,
            sam=sam,
            fwd_reads=fwd_reads,
            rev_reads=rev_reads,
            se_reads=se_reads,
            cpus=cpus,
        ):
            bowtie.build(fasta, db)
            bowtie.align(db, sam, fwd_reads, rev_reads, se_reads, cpus=cpus)

        # Setup of coverage calculation sequence depending on file(s) provided
        calculation_sequence = {
            "bed_exists": [parse_bed],
            "bam_exists": [make_bed, parse_bed],
            "sam_exists": [sort_samfile, make_bed, parse_bed],
            "full": [align_reads, sort_samfile, make_bed, parse_bed],
        }
        # Now we need to determine which point to start the calculation...
        for fp, argname in zip([bed, bam, sam], ["bed", "bam", "sam"]):
            step = "full"
            if os.path.exists(fp):
                step = f"{argname}_exists"
                break

        if (not fwd_reads or not rev_reads) and step == "full":
            raise ValueError(
                f"fwd_reads and rev_reads are required if no other alignments are specified!"
            )
        logger.debug(f"starting coverage calculation sequence from {step}")
        for calculation in calculation_sequence[step]:
            logger.debug(f"running {calculation.__name__}")
            if calculation.__name__ == "parse_bed":
                return calculation()
            calculation()


def main():
    import argparse
    import logging as logger
    import multiprocessing as mp

    logger.basicConfig(
        format="%(asctime)s : %(name)s : %(levelname)s : %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
    Construct contig coverage table given an input `assembly` and provided files.

    Provided files may include one from the list below:
    1. `fwd_reads` and/or `rev_reads` and/or `se_reads`
    2. `sam` - alignment of `assembly` and `reads` in SAM format
    3. `bam` - alignment of `assembly` and `reads` in BAM format
    4. `bed` - alignment of `assembly` and `reads` in BED format
    """,
    )
    parser.add_argument(
        "-f", "--assembly", help="</path/to/metagenome.fasta>", required=True
    )
    parser.add_argument(
        "-1", "--fwd-reads", help="</path/to/forwards-reads.fastq>", nargs="*"
    )
    parser.add_argument(
        "-2", "--rev-reads", help="</path/to/reverse-reads.fastq>", nargs="*"
    )
    parser.add_argument(
        "-U", "--se-reads", help="</path/to/single-end-reads.fastq>", nargs="*"
    )
    parser.add_argument("--sam", help="</path/to/alignments.sam>")
    parser.add_argument("--bam", help="</path/to/alignments.bam>")
    parser.add_argument("--bed", help="</path/to/alignments.bed>")
    parser.add_argument(
        "--cpus",
        help=f"Num processors to use. (default: {mp.cpu_count()})",
        default=mp.cpu_count(),
        type=int,
    )
    parser.add_argument(
        "--from-spades",
        help="Extract k-mer coverages from contig IDs. (Input assembly is output from SPAdes)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--out", help="Path to write a table of coverages", required=True
    )
    args = parser.parse_args()

    if args.from_spades:
        fh = (
            gzip.open(args.assembly, "rt")
            if args.assembly.endswith(".gz")
            else open(args.assembly)
        )
        records = [rec for rec in SeqIO.parse(fh, "fasta")]
        fh.close()
        coverages = from_spades_names(records)
        logger.info(
            f"{coverages.index.nunique():,} contig coverages retrieved from {args.assembly}"
        )
        coverages.to_csv(args.out, sep="\t", index=True, header=True)
        logger.info(f"written: {args.out}")
        return

    get(
        fasta=args.assembly,
        fwd_reads=args.fwd_reads,
        rev_reads=args.rev_reads,
        sam=args.sam,
        bam=args.bam,
        bed=args.bed,
        cpus=args.cpus,
        out=args.out,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        logger.exception(err)
        import sys

        print("Coverage calculation failed. check the logs for details")
        sys.exit(1)
