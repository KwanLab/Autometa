#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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
Class and functions related to running diamond on metagenome sequences
"""


import logging
import os
import subprocess

import multiprocessing as mp

from tqdm import tqdm

from autometa.common.utilities import file_length

logger = logging.getLogger(__name__)


def makedatabase(fasta: str, database: str, cpus: int = mp.cpu_count()) -> str:
    """
    Creates a database against which the query sequence would be blasted

    Parameters
    ----------
    fasta : str
        Path to fasta file whose database needs to be made
        e.g. '<path/to/fasta/file>'
    database : str
        Path to the output diamond formatted database file
        e.g. '<path/to/database/file>'
    cpus : int, optional
        Number of processors to be used. By default uses all the processors of the system

    Returns
    -------
    str
        Path to diamond formatted database

    Raises
    ------
    subprocess.CalledProcessError
        Failed to create diamond formatted database
    """
    cmd = ["diamond", "makedb", "--in", fasta, "--db", database, "-p", str(cpus)]
    logger.debug(" ".join(cmd))
    subprocess.run(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )
    return database


def blast(
    fasta: str,
    database: str,
    outfpath: str,
    blast_type: str = "blastp",
    evalue: float = float("1e-5"),
    maxtargetseqs: int = 200,
    cpus: int = mp.cpu_count(),
    tmpdir: str = None,
    force: bool = False,
    verbose: bool = False,
) -> str:
    """
    Performs diamond blastp search using query sequence against diamond formatted database

    Parameters
    ----------
    fasta : str
        Path to fasta file having the query sequences. Should be amino acid sequences in case of BLASTP
        and nucleotide sequences in case of BLASTX
    database : str
        Path to diamond formatted database
    outfpath : str
        Path to output file
    blast_type : str, optional
        blastp to align protein query sequences against a protein reference database,
        blastx to align translated DNA query sequences against a protein reference database, by default 'blastp'
    evalue : float, optional
        cutoff e-value to count hit as significant, by default float('1e-5')
    maxtargetseqs : int, optional
        max number of target sequences to retrieve per query by diamond, by default 200
    cpus : int, optional
        Number of processors to be used, by default uses all the processors of the system
    tmpdir : str, optional
        Path to temporary directory. By default, same as the output directory
    force : bool, optional
        overwrite existing diamond results, by default False
    verbose : bool, optional
        log progress to terminal, by default False

    Returns
    -------
    str
        Path to BLAST results

    Raises
    ------
    FileNotFoundError
        `fasta` file does not exist
    ValueError
        provided `blast_type` is not 'blastp' or 'blastx'
    subprocess.CalledProcessError
        Failed to run blast
    """
    if not os.path.exists(fasta):
        raise FileNotFoundError(fasta)
    if os.path.exists(outfpath) and not force:
        if os.path.getsize(outfpath):
            if verbose:
                logger.warning(f"FileExistsError: {outfpath}. To overwrite use --force")
            return outfpath
    blast_type = blast_type.lower()
    if blast_type not in ["blastp", "blastx"]:
        raise ValueError(f"blast_type must be blastp or blastx. Given: {blast_type}")
    if verbose:
        logger.debug(f"diamond {blast_type} {fasta} against {database}")
    cmd = [
        "diamond",
        blast_type,
        "--query",
        fasta,
        "--db",
        database,
        "--evalue",
        evalue,
        "--max-target-seqs",
        maxtargetseqs,
        "--threads",
        cpus,
        "--outfmt",
        "6",
        "--out",
        outfpath,
    ]
    if tmpdir:
        cmd.extend(["--tmpdir", tmpdir])
    # this is done as when cmd is a list each element should be a string
    cmd = [str(c) for c in cmd]
    if verbose:
        logger.debug(f'RunningDiamond: {" ".join(cmd)}')
    subprocess.run(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )
    return outfpath


def parse(results: str, bitscore_filter: float = 0.9, verbose: bool = False) -> dict:
    """
    Retrieve diamond results from output table

    Parameters
    ----------
    results : str
        Path to BLASTP output file in outfmt9
    bitscore_filter : 0 < float <= 1, optional
        Bitscore filter applied to each sseqid, by default 0.9
        Used to determine whether the bitscore is above a threshold value.
        For example, if it is 0.9 then only bitscores >= 0.9 * the top bitscore are accepted
    verbose : bool, optional
        log progress to terminal, by default False

    Returns
    -------
    dict
        {qseqid: {sseqid, sseqid, ...}, ...}

    Raises
    -------
    FileNotFoundError
        diamond results table does not exist
    ValueError
        bitscore_filter value is not a float or not in range of 0 to 1
    """
    disable = False if verbose else True
    # boolean toggle --> keeping above vs. below because I think this is more readable.
    # disable = not verbose
    if verbose:
        logger.debug(f"Parsing accessions from {os.path.basename(results)}")
    if not os.path.exists(results):
        raise FileNotFoundError(results)
    try:
        float(bitscore_filter)
    except ValueError:
        raise ValueError(
            f"bitscore_filter must be a float! Input: {bitscore_filter} Type: {type(bitscore_filter)}"
        )
    in_range = 0.0 < bitscore_filter <= 1.0
    if not in_range:
        raise ValueError(f"bitscore_filter not in range(0,1)! Input: {bitscore_filter}")
    hits = {}
    n_lines = file_length(results) if verbose else None
    topbitscore = float("-inf")
    with open(results) as fh:
        for line in tqdm(
            fh, disable=disable, total=n_lines, desc="Parsing Accessions", leave=False
        ):
            llist = line.rstrip().split("\t")
            qseqid = llist[0]
            sseqid = llist[1]
            bitscore = float(llist[11])
            # Reassign the topbitscore if this is a new qseqid from the BLAST table.
            if qseqid not in hits:
                hits.update({qseqid: set([sseqid])})
                topbitscore = bitscore
                continue
            if bitscore >= bitscore_filter * topbitscore:
                hits[qseqid].add(sseqid)
    return hits


if __name__ == "__main__":
    # fmt: off
    import sys
    print("Diamond contains utility functions and should not be run directly!")
    sys.exit(0)
