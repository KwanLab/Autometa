#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Class and functions related to running diamond on metagenome sequences
"""


import gzip
import logging
import os
import subprocess

from itertools import chain
from tqdm import tqdm

from autometa.common.utilities import make_pickle,file_length

logger = logging.getLogger(__name__)


class DiamondResult:
    """docstring for DiamondResult.

    Some operator overloading here...
    For other examples of this see:
        https://www.geeksforgeeks.org/operator-overloading-in-python/

    """
    def __init__(self, qseqid, sseqid, pident, length, mismatch, gapopen,
        qstart, qend, sstart, send, evalue, bitscore):
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.pident = float(pident)
        self.length = int(length)
        self.mismatch = int(mismatch)
        self.gapopen = int(gapopen)
        self.qstart = int(qstart)
        self.qend = int(qend)
        self.sstart = int(sstart)
        self.send = int(send)
        self.evalue = float(evalue)
        self.bitscore = float(bitscore)
        self.sseqids = {self.sseqid:{
            'pident':self.pident,
            'length':self.length,
            'mismatch':self.mismatch,
            'gapopen':self.gapopen,
            'qstart':self.qstart,
            'qend':self.qend,
            'sstart':self.sstart,
            'send':self.send,
            'evalue':self.evalue,
            'bitscore':self.bitscore,
        }}

    # def __repr__(self):
    #     return str(self)

    def __str__(self):
        return f'''{self.qseqid}
Num. sseqids: {len(self.sseqids)}
Top sseqid: {self.get_top_hit()}
'''

    def __eq__(self, other_hit):
        if self.qseqid == other_hit.qseqid:
            return True
        else:
            return False

    def __add__(self, other_hit):
        assert self == other_hit, f'qseqids do not match! {self} & {other_hit}'
        self.sseqids.update(other_hit.sseqids)
        return self

    def __sub__(self, other_hit):
        assert self == other_hit, f'qseqids do not match! {self} & {other_hit}'
        self.sseqids.pop(other_hit.sseqid)
        return self

    def get_top_hit(self):
        top_bitscore = float('-Inf')
        for sseqid,attrs in self.sseqids.items():
            if attrs.get('bitscore') >= top_bitscore:
                top_bitscore = attrs.get('bitscore')
                top_hit = sseqid
        return top_hit

def makedatabase(fasta, database, nproc=1):
    cmd = f'diamond makedb --in {fasta} --db {database} -p {nproc}'
    logger.debug(f'{cmd}')
    with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
        retcode = subprocess.call(cmd, stdout=stdout, stderr=stderr, shell=True)
    if retcode:
        raise OSError(f'DiamondFailed:\nArgs:{proc.args}\nReturnCode:{proc.returncode}')
    return database

def blast(fasta, database, outfpath, blast_type='blastp', evalue=float('1e-5'),
    maxtargetseqs=200, cpus=0, tmpdir=os.curdir, force=False, verbose=False):
    """Performs diamond blastp search using fasta against diamond formatted database

    Parameters
    ----------
    fasta : str
        </path/to/fasta/file>. May be amino acid or nucleotide sequences
    database : str
        </path/to/diamond/formatted/database>.
    outfpath : str
        </path/to/output/file>.
    blast_type : str
        blastp if fasta consists of amino acids or blastx if nucleotides
    evalue : float
        cutoff e-value to count hit as significant (the default is float('1e-5')).
    maxtargetseqs : int
        max number of target sequences to retrieve per query by diamond (the default is 200).
    cpus : int
        number of cpus to use (the default is 1).
    tmpdir : type
        </path/to/temporary/directory> (the default is os.curdir).
    force : boolean
        overwrite existing diamond results `force` (the default is False).
    verbose : boolean
        log progress to terminal `verbose` (the default is False).

    Returns
    -------
    str
        `outfpath`

    Raises
    -------
    FileNotFoundError
        `fasta` file does not exist
    ValueError
        provided `blast_type` is not 'blastp' or 'blastx'
    OSError
        Diamond execution failed
    """
    if not os.path.exists(fasta):
        raise FileNotFoundError(fasta)
    if os.path.exists(outfpath) and not force:
        empty = not os.stat(outfpath).st_size
        if not empty:
            if verbose:
                logger.warning(f'FileExistsError: {outfpath}. To overwrite use --force')
            return outfpath
    blast_type = blast_type.lower()
    if blast_type not in ['blastp','blastx']:
        raise ValueError(f'blast_type must be blastp or blastx. Given: {blast_type}')
    if verbose:
        logger.debug(f'Diamond{blast_type.title()} {fasta} against {database}')
    cmd = [
        'diamond',
        'blastp',
        '--query',
        fasta,
        '--db',
        database,
        '--evalue',
        evalue,
        '--max-target-seqs',
        maxtargetseqs,
        '--threads',
        cpus,
        '--outfmt',
        '6',
        '--out',
        outfpath,
        '--tmpdir',
        tmpdir]
    cmd = [str(c) for c in cmd]
    if verbose:
        logger.debug(f'RunningDiamond: {" ".join(cmd)}')
    with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
        proc = subprocess.run(cmd, stdout=stdout, stderr=stderr)
    if proc.returncode:
        raise OSError(f'DiamondFailed:\nArgs:{proc.args}\nReturnCode:{proc.returncode}')
    return outfpath

def parse(results, top_pct=0.9, verbose=False):
    """Retrieve diamond results from output table

    Parameters
    ----------
    results : str
        </path/to/blastp/outfmt6/output/file>
    top_pct : 0 < float <= 1
        bitscore filter applied to each qseqid (the default is 0.9).

    Returns
    -------
    dict
        {qseqid: DiamondResult, ...}

    Raises
    -------
    FileNotFoundError
        diamond results table does not exist
    ValueError
        top_pct value is not a float or not in range of 0 to 1
    """
    disable = False if verbose else True
    # boolean toggle --> keeping above vs. below because I think this is more readable.
    # disable = not verbose
    if verbose:
        logger.debug(f'Parsing accessions from {os.path.basename(results)}')
    if not os.path.exists(results):
        raise FileNotFoundError(results)
    try:
        float(top_pct)
    except ValueError as err:
        raise ValueError(f'top_pct must be a float! Input: {top_pct} Type: {type(top_pct)}')
    in_range = 0.0 < top_pct <= 1.0
    if not in_range:
        raise ValueError(f'top_pct not in range(0,1)! Input: {top_pct}')
    hits = {}
    temp = set()
    n_lines = file_length(results) if verbose else None
    with open(results) as fh:
        for line in tqdm(fh, disable=disable, total=n_lines, desc='Parsing Accessions',leave=False):
            llist = line.rstrip().split('\t')
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
            evalue = llist[10]
            bitscore = llist[11]
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
                bitscore=bitscore)
            if hit.qseqid not in temp:
                hits.update({hit.qseqid:hit})
                topbitscore = hit.bitscore
                temp = set([hit.qseqid])
                continue
            if hit.bitscore >= top_pct * topbitscore:
                hits[hit.qseqid] += hit
    return hits

def add_taxids(hits, database, verbose=True):
    """Translates accessions to taxid translations from prot.accession2taxid.gz

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
        </path/to/prot.accession2taxid.gz>

    Returns
    -------
    dict
        hits - {qseqid: DiamondResult, ...} sseqids dicts updated with 'taxid' key

    Raises
    -------
    FileNotFoundError
        prot.accession2taxid.gz database is required for translation taxid
    DatabasesOutOfDateError
        prot.accession2taxid.gz database and nr.dmnd are out of sync resulting
        in accessions that are no longer available (either due to being
        suppressed, deprecated or removed by NCBI). This must be resolved by
        updating both nr and prot.accession2taxid.gz and re-running diamond on
        the new nr.dmnd database. Alternatively, can try to find the exceptions in merged.dmp


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
    accessions = set(chain.from_iterable([hit.sseqids.keys() for qseqid,hit in hits.items()]))
    fh = gzip.open(database) if database.endswith('.gz') else open(database)
    __ = fh.readline()
    if verbose:
        logger.debug(f'Searching for {len(accessions):,} accessions in {os.path.basename(database)}. This may take a while...')
    is_gzipped = True if database.endswith('.gz') else False
    n_lines = file_length(database) if verbose else None
    desc = f'Parsing {os.path.basename(database)}'
    acc2taxids = {}
    for line in tqdm(fh, disable=disable, desc=desc, total=n_lines, leave=False):
        line = line.decode() if is_gzipped else line
        acc_num, acc_ver, taxid, _ = line.split('\t')
        taxid = int(taxid)
        if acc_num in accessions:
            acc2taxids.update({acc_num:taxid})
        if acc_ver in accessions:
            acc2taxids.update({acc_ver:taxid})
    fh.close()
    n_qseqids = len(hits)
    desc = f'Translating {n_qseqids} qseqids\' accessions to taxids'
    for qseqid,hit in tqdm(hits.items(), disable=disable, total=n_qseqids, desc=desc, leave=False):
        for sseqid in hit.sseqids:
            taxid = acc2taxids.get(sseqid)
            # if taxid is None:
            #     raise DatabasesOutOfDateError(f'{sseqid} deprecated/suppressed/removed')
            hit.sseqids[sseqid].update({'taxid':taxid})
    return hits

def main(args):
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
        verbose=args.verbose)
    hits = parse(results=result, top_pct=args.top_pct, verbose=args.verbose)
    hits = add_taxids(
        hits=hits,
        database=args.acc2taxids,
        verbose=args.verbose)
    fname,__ = os.path.splitext(os.path.basename(args.outfile))
    dirpath = os.path.dirname(os.path.realpath(args.outfile))
    hits_fname = '.'.join([fname, 'pkl.gz'])
    hits_fpath = os.path.join(dirpath, hits_fname)
    pickled_fpath = make_pickle(obj=hits, outfpath=hits_fpath)
    logger.debug(f'{len(hits):,} diamond hits serialized to {pickled_fpath}')

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Retrieves blastp hits with provided input assembly')
    parser.add_argument('fasta', help='</path/to/faa/file>')
    parser.add_argument('database', help='</path/to/diamond/formatted/database>')
    parser.add_argument('acc2taxids', help='</path/to/ncbi/prot.accession2taxid.gz>')
    parser.add_argument('outfile', help='</path/to/diamond/output/file>')
    parser.add_argument('blast_type', help='[blastp]: A.A -> A.A. [blastx]: Nucl. -> A.A.',
        default='blastp',choices=['blastp','blastx'])
    parser.add_argument('--evalue', help='diamond evalue threshold', default=float('1e-5'))
    parser.add_argument('--maxtargetseqs',
        help='max target sequences to retrieve per query', default=200, type=int)
    parser.add_argument('--cpus', help='num cpus to use', default=0, type=int)
    parser.add_argument('--tmpdir', help='</path/to/tmp/directory>', default=os.curdir)
    parser.add_argument('--top-pct', help='top percentage of hits to retrieve', default=0.9)
    parser.add_argument('--force', help='force overwrite of diamond output table',
        action='store_true')
    parser.add_argument('--verbose', help='add verbosity', action='store_true')
    args = parser.parse_args()
    main(args)
