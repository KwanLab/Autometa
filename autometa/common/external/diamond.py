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
Class and functions related to running diamond on metagenome sequences
"""


import gzip
import logging
import os
import subprocess

import multiprocessing as mp

from itertools import chain
from tqdm import tqdm

from autometa.common.utilities import make_pickle, file_length
from autometa.common.external import prodigal
from autometa.common.exceptions import DatabaseOutOfSyncError

logger = logging.getLogger(__name__)


class DiamondResult:
    """
    DiamondResult class to instantiate a DiamondResult object for each qseqid.
    Allows for better handling of subject sequence IDs that hit to a qseqid.

    Includes methods used to modify the DiamondResult object (add or remove sseqid from sseqids dictionary),
    check if two DiamondResult objects have the same qseqid and return a user friendly and unambiguous output
    from str() and repr() respectively. Also used to return the sseqid with the highest bitscore amongst all
    the subject sequences that hit a query. These methods come in handy when retrieving diamond results from output table.

    Attributes
    ----------
    sseqids : dict
        All the subject sequences that hit to the query sequence
        {sseqid:parameters, sseqid:parameters, ...}
    qseqid: str
        result query sequence ID

    """

    def __init__(
        self,
        qseqid,
        sseqid,
        pident,
        length,
        mismatch,
        gapopen,
        qstart,
        qend,
        sstart,
        send,
        evalue,
        bitscore,
    ):
        """
        Instantiates the DiamondResult class

        Parameters
        ----------
        qseqid : str
            query sequence ID
        sseqid : str
            subject sequence ID
        pident : float
            Percentage of identical matches.
        length : int
            Alignment length.
        mismatch : int
            Number of mismatches.
        gapopen : int
            Number of gap openings.
        qstart : int
            Start of alignment in query.
        qend : int
            End of alignment in query.
        sstart : int
            Start of alignment in subject.
        send : int
            End of alignment in subject sequence.
        evalue : float
            Expect value.
        bitscore : float
            Bitscore.
        """
        self.qseqid = qseqid
        self.sseqids = {
            sseqid: {
                "pident": float(pident),
                "length": int(length),
                "mismatch": int(mismatch),
                "gapopen": int(gapopen),
                "qstart": int(qstart),
                "qend": int(qend),
                "sstart": int(sstart),
                "send": int(send),
                "evalue": float(evalue),
                "bitscore": float(bitscore),
            }
        }

    def __repr__(self):
        """
        Operator overloading to return the representation of the class object

        Returns
        -------
        str
            Class name followed by query sequence ID
        """
        return f"Class: {self.__class__.__name__}, Query seqID: {self.qseqid}"

    def __str__(self):
        """
        Operator overloading to return the string representation of the class objects

        Returns
        -------
        str
            String representation of query sequence ID, followed by total number of hits and finally the highest
            bit score of all the hits
        """
        return f"{self.qseqid}; {len(self.sseqids)} sseqids; top hit by bitscore: {self.get_top_hit()}"

    def __eq__(self, other_hit):
        """
        Operator overloading to compare two objects of the DiamondResult class

        Parameters
        ----------
        other_hit : DiamondResult object
            Other DiamondResult object to compare with

        Returns
        -------
        Boolean
            True if qseqid corresponding to both DiamondResult objects are equal, else False
        """
        if self.qseqid == other_hit.qseqid:
            return True
        else:
            return False

    def __add__(self, other_hit):
        """
        Operator overloading to update (add) the sseqids dictionary with the other sseqid hit dictionary.
        The addition will only be successful if both the DiamondResult objects have the same qseqid

        Parameters
        ----------
        other_hit : DiamondResult object
            Other DiamondResult object to add

        Returns
        -------
        DiamondResult object
            DiamondResult object whose sseqids dict has been updated (added) with another sseqid

        Raises
        ------
        AssertionError
            Query sequences are not equal
        """
        assert (
            self == other_hit
        ), f"qseqids do not match! {self.qseqid} & {other_hit.qseqid}"
        self.sseqids.update(other_hit.sseqids)
        return self

    def __sub__(self, other_hit):
        """
        Operator overloading to remove (subtract) the other sseqid hits dictionary from the sseqids dictionary.
        The subtraction will only be successful if both the DiamondResult objects have the same qseqid

        Parameters
        ----------
        other_hit : Dict
            Other DiamondResult object to subtract

        Returns
        -------
        DiamondResult object
            DiamondResult object where a sseqid has been removed (subtracted) from the sseqids dict

        Raises
        ------
        AssertionError
            Query sequences are not equal
        """
        assert (
            self == other_hit
        ), f"qseqids do not match! {self.qseqid} & {other_hit.qseqid}"
        for sseqid in other_hit.sseqids:
            try:
                self.sseqids.pop(sseqid)
            except KeyError:
                raise KeyError(
                    f"Given sseqid: {sseqid} is absent from the corresponding DiamondResult"
                )
        return self

    def get_top_hit(self):
        """
        Returns the subject sequence ID with the highest bitscore amongst all the subject sequences that hit a query
        """
        top_bitscore = float("-Inf")
        top_hit = None
        for sseqid, attrs in self.sseqids.items():
            if attrs.get("bitscore", float("-Inf")) >= top_bitscore:
                top_bitscore = attrs.get("bitscore")
                top_hit = sseqid
        return top_hit


def makedatabase(fasta, database, cpus=mp.cpu_count()):
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
    fasta,
    database,
    outfpath,
    blast_type="blastp",
    evalue=float("1e-5"),
    maxtargetseqs=200,
    cpus=mp.cpu_count(),
    tmpdir=None,
    force=False,
    verbose=False,
):
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


def parse(results, bitscore_filter=0.9, verbose=False):
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
        {qseqid: DiamondResult, ...}

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
    temp = set()
    n_lines = file_length(results) if verbose else None
    with open(results) as fh:
        for line in tqdm(
            fh, disable=disable, total=n_lines, desc="Parsing Accessions", leave=False
        ):
            llist = line.rstrip().split("\t")
            qseqid = llist[0]
            sseqid = llist[1]
            pident = llist[2]
            length = llist[3]
            mismatch = llist[4]
            gapopen = llist[5]
            qstart = llist[6]
            qend = llist[7]
            sstart = llist[8]
            send = llist[9]
            evalue = float(llist[10])
            bitscore = float(llist[11])
            hit = DiamondResult(
                qseqid=qseqid,
                sseqid=sseqid,
                pident=pident,
                length=length,
                mismatch=mismatch,
                gapopen=gapopen,
                qstart=qstart,
                qend=qend,
                sstart=sstart,
                send=send,
                evalue=evalue,
                bitscore=bitscore,
            )
            if hit.qseqid not in temp:
                hits.update({hit.qseqid: hit})
                topbitscore = bitscore
                temp = set([hit.qseqid])
                continue
            if bitscore >= bitscore_filter * topbitscore:
                hits[hit.qseqid] += hit
    return hits


def add_taxids(hits, database, verbose=True):
    """
    Translates accessions to taxids from prot.accession2taxid.gz. If an accession number is no
    longer available in prot.accesssion2taxid.gz (either due to being suppressed, deprecated or
    removed by NCBI), then None is returned as the taxid for the corresponsing sseqid.


    # TODO: Should maybe write a wrapper for this to run on all of the NCBI
    databases listed below... Maybe this will help account for instances where the
    accession is not found due to suppression,deprecation,removal,etc...

    DiamondResults contain sseqids with non-redundant protein sequence database
    entries from:
    - GenPept
    - Swissprot
    - PIR
    - PDF
    - PDB
    - RefSeq

    Parameters
    ----------
    hits : dict
        {qseqid: DiamondResult, ...}
    database : str
        Path to prot.accession2taxid.gz database
    verbose : bool, optional
        log progress to terminal, by default False

    Returns
    -------
    dict
        hits - {qseqid: DiamondResult, ...} sseqids dicts updated with 'taxid' key

    Raises
    -------
    FileNotFoundError
        prot.accession2taxid.gz database is required for translation taxid

    # TODO: Replace file_length func for database file.
    (in this case 808,717,857 lines takes ~15 minutes simply to read each line...)
    This causes unnecessary slowdown. On download of the respective database. Line
    counts can be saved in config and replaced here for the file_length func. This
    will allow for only one read-through of the file, instead of two.
    """
    disable = not verbose
    if not os.path.exists(database):
        raise FileNotFoundError(database)
    # OPTIMIZE: prot.accession2taxid.gz database is almost 1 billion lines...
    accessions = set(
        chain.from_iterable([hit.sseqids.keys() for qseqid, hit in hits.items()])
    )
    # "rt" open the database in text mode instead of binary. Now it can be handled like a text file
    fh = gzip.open(database, "rt") if database.endswith(".gz") else open(database)
    __ = fh.readline()  # remove the first line as it just gives the description
    if verbose:
        logger.debug(
            f"Searching for {len(accessions):,} accessions in {os.path.basename(database)}. This may take a while..."
        )
    n_lines = file_length(database, approximate=True) if verbose else None
    desc = f"Parsing {os.path.basename(database)}"
    acc2taxids = {}
    for line in tqdm(fh, disable=disable, desc=desc, total=n_lines, leave=False):
        acc_num, acc_ver, taxid, _ = line.split("\t")
        taxid = int(taxid)
        if acc_num in accessions:
            acc2taxids.update({acc_num: taxid})
        if acc_ver in accessions:
            acc2taxids.update({acc_ver: taxid})
    fh.close()
    n_qseqids = len(hits)
    desc = f"Translating {n_qseqids} qseqids' accessions to taxids"
    for qseqid, hit in tqdm(
        hits.items(), disable=disable, total=n_qseqids, desc=desc, leave=False
    ):
        for sseqid in hit.sseqids:
            taxid = acc2taxids.get(sseqid)
            if not taxid:
                raise DatabaseOutOfSyncError(
                    f"{sseqid} is either a new taxid not present in prot.accession2taxid or is deprecated/suppressed/removed"
                )
            hit.sseqids[sseqid].update({"taxid": taxid})
    return hits


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="""
    Retrieves blastp hits with provided input assembly
    """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("fasta", help="Path to fasta file having the query sequences")
    parser.add_argument("database", help="Path to diamond formatted database")
    parser.add_argument("acc2taxids", help="Path to prot.accession2taxid.gz database")
    parser.add_argument("outfile", help="Path to output file")
    parser.add_argument(
        "blast_type",
        help="[blastp]: A.A -> A.A. [blastx]: Nucl. -> A.A.",
        default="blastp",
        choices=["blastp", "blastx"],
    )
    parser.add_argument(
        "--evalue", help="diamond evalue threshold", default=float("1e-5")
    )
    parser.add_argument(
        "--maxtargetseqs",
        help="max target sequences to retrieve per query",
        default=200,
        type=int,
    )
    parser.add_argument(
        "--cpus", help="number of processors to use", default=mp.cpu_count(), type=int
    )
    parser.add_argument(
        "--tmpdir",
        help="Path to directory which will be used for temporary storage by diamond",
    )
    parser.add_argument(
        "--bitscore-filter",
        help="hits with bitscore greater than bitscr-filter*top bitscore will be kept",
        default=0.9,
    )
    parser.add_argument(
        "--force", help="force overwrite of diamond output table", action="store_true"
    )
    parser.add_argument("--verbose", help="add verbosity", action="store_true")
    args = parser.parse_args()
    result = blast(
        fasta=args.fasta,
        database=args.database,
        outfpath=args.outfile,
        blast_type=args.blast_type,
        evalue=args.evalue,
        maxtargetseqs=args.maxtargetseqs,
        cpus=args.cpus,
        tmpdir=args.tmpdir,
        force=args.force,
        verbose=args.verbose,
    )
    hits = parse(
        results=result, bitscore_filter=args.bitscore_filter, verbose=args.verbose
    )
    hits = add_taxids(hits=hits, database=args.acc2taxids, verbose=args.verbose)
    fname, __ = os.path.splitext(os.path.basename(args.outfile))
    dirpath = os.path.dirname(os.path.realpath(args.outfile))
    hits_fname = ".".join([fname, "pkl.gz"])
    hits_fpath = os.path.join(dirpath, hits_fname)
    pickled_fpath = make_pickle(obj=hits, outfpath=hits_fpath)
    logger.debug(f"{len(hits):,} diamond hits serialized to {pickled_fpath}")


if __name__ == "__main__":
    main()
