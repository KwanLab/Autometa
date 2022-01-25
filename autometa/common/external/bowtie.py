#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
Script containing wrapper functions for bowtie2.
"""


import logging
import os
import subprocess
from typing import List


logger = logging.getLogger(__name__)


def run(cmd: str) -> bool:
    """Run `cmd` via subprocess.

    Parameters
    ----------
    cmd : str
        Executable input str

    Returns
    -------
    bool
        True if no returncode from subprocess.call else False

    """
    logger.debug(f"run: {cmd}")
    with open(os.devnull, "w") as stdout, open(os.devnull, "w") as stderr:
        retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
    if retcode:
        logger.warning(f"Args:{cmd} ReturnCode:{retcode}")
        return False
    return True


def build(assembly: str, out: str) -> str:
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
    ChildProcessError
        bowtie2-build failed
    """
    cmd = f"bowtie2-build {assembly} {out}"
    success = run(cmd)
    if not success:
        raise ChildProcessError(f"{cmd} failed. {out} not written")
    return out


def align(
    db: str,
    sam: str,
    fwd_reads: List[str] = None,
    rev_reads: List[str] = None,
    se_reads: List[str] = None,
    cpus: int = 0,
    **kwargs,
) -> str:
    """Align reads to bowtie2-index `db` (at least one `*_reads` argument is required).

    Parameters
    ----------
    db : str
        </path/to/prefix/bowtie2/database>. I.e. `db`.{#}.bt2
    sam : str
        </path/to/out.sam>
    fwd_reads : list, optional
        [</path/to/forward_reads.fastq>, ...]
    rev_reads : list, optional
        [</path/to/reverse_reads.fastq>, ...]
    se_reads : list, optional
        [</path/to/single_end_reads.fastq>, ...]
    cpus : int, optional
        Num. processors to use (the default is 0).
    **kwargs : dict, optional
        Additional optional args to supply to bowtie2. Must be in format:
        key = flag
        value = flag-value

    Returns
    -------
    str
        </path/to/out.sam>

    Raises
    -------
    ChildProcessError
        bowtie2 failed
    """
    exe = f"bowtie2 -x {db}"
    flags = "-q --phred33 --very-sensitive --no-unal"
    sam_out = f"-S {sam}"
    params = [exe, flags, sam_out]
    if not isinstance(cpus, int) or cpus < 0:
        raise ValueError(f"cpus must be an integer greater than 0. given: {cpus}")
    # cpus==0 will skip adding -p/--threads flag
    if cpus:
        params.append(f"-p {cpus}")
    reads_provided = False
    for flag, reads in zip(["-1", "-2", "-U"], [fwd_reads, rev_reads, se_reads]):
        if reads:
            reads_provided = True
            if isinstance(reads, list):
                reads_to_str = ",".join([str(read) for read in reads])
                params.append(f"{flag} {reads_to_str}")
            else:
                params.append(f"{flag} {reads}")
    if not reads_provided:
        raise ValueError(f"At least one fastq file is required!")
    if kwargs:
        params += [f"{flag} {value}" for flag, value in kwargs.items()]
    cmd = " ".join(params)
    success = run(cmd)
    if not success:
        raise ChildProcessError(f"{cmd} failed. {sam} not written")
    return sam


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Align provided reads to metagenome `assembly` and write alignments to `sam`."
        "NOTE: At least one reads file is required."
    )
    parser.add_argument("assembly", help="</path/to/assembly.fasta>")
    parser.add_argument(
        "database",
        help="</path/to/alignment.database>. Will construct database at provided path if not found.",
    )
    parser.add_argument("sam", help="</path/to/alignment.sam>")
    parser.add_argument(
        "-1", "--fwd-reads", help="</path/to/forward-reads.fastq>", nargs="*"
    )
    parser.add_argument(
        "-2", "--rev-reads", help="</path/to/reverse-reads.fastq>", nargs="*"
    )
    parser.add_argument(
        "-U", "--se-reads", help="</path/to/single-end-reads.fastq>", nargs="*"
    )
    parser.add_argument("--cpus", help="Num processors to use.", default=1, type=int)
    args = parser.parse_args()

    build(args.assembly, args.database)
    align(
        database=args.database,
        sam=args.sam,
        fwd_reads=args.fwd_reads,
        rev_reads=args.rev_reads,
        se_reads=args.se_reads,
        cpus=args.cpus,
        kwargs=args.kwargs,
    )


if __name__ == "__main__":
    main()
