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

Autometa Bin Class
"""


import logging
import os

import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqUtils
from functools import lru_cache

from autometa.common import kmers
from autometa.common.utilities import timeit
from autometa.common.external import prodigal


logger = logging.getLogger(__name__)


class MetaBin:
    """Autometa MetaBin class to manipulate/annotate contigs from a metagenome.

    Parameters
    ----------
    assembly : str
        </path/to/metagenome.fasta>
    contig_ids : list
        List of contig ids to manipulate/annotate (must be contained in metagenome).
    seqrecords : list
        List of seqrecords to manipulate/annotate (must be contained in metagenome).
    outdir : str, optional
        </path/to/output/directory> (Default is the directory storing the `assembly`).

    Attributes
    ----------
    basename : str
        base name of `assembly`
    root : str
        root name of `assembly` (Will remove common extension like '.fasta')
    nucls_fname : str
        File name of contigs nucleotide ORFs
    prots_fname : str
        File name of contigs amino-acid ORFs
    nucl_orfs_fpath : str
        </path/to/nucleotide/`outdir`/`nucls_fname`>
    prot_orfs_fpath : str
        </path/to/amino/acid/`outdir`/`prots_fname`>
    nseqs : int
        Number of contigs in MetaBin

    Raises
    ------
    ValueError
        One of `contig_ids` or `seqrecords` must be provided
    ValueError
        `contig_ids` do not match `seqrecords`
    """

    def __init__(self, assembly, contig_ids=[], seqrecords=[], outdir=None):
        if not contig_ids and not seqrecords:
            raise ValueError("One of `contig_ids` or `seqrecords` must be provided!")
        if contig_ids and seqrecords:
            contig_id_diff = set(contig_ids).difference({rec.id for rec in seqrecords})
            if contig_id_diff:
                raise ValueError(
                    f"`contig_ids` do not match `seqrecords`. Extra `contig_ids`: {contig_id_diff}"
                )
            seqrecord_diff = {rec.id for rec in seqrecords}.difference(contig_ids)
            if seqrecord_diff:
                raise ValueError(
                    f"`seqrecords` do not match `contig_ids`. Extra `seqrecords`: {seqrecord_diff}"
                )
        self.assembly = os.path.realpath(assembly)
        self.basename = os.path.basename(self.assembly)
        self.root = os.path.splitext(self.basename)[0]
        self.outdir = (
            os.path.realpath(outdir) if outdir else os.path.dirname(self.assembly)
        )
        nucls_ext = "orfs.fna"
        prots_ext = "orfs.faa"
        self.nucls_fname = ".".join([self.root, nucls_ext])
        self.prots_fname = ".".join([self.root, prots_ext])
        self.nucl_orfs_fpath = os.path.join(self.outdir, self.nucls_fname)
        self.prot_orfs_fpath = os.path.join(self.outdir, self.prots_fname)
        if contig_ids and not seqrecords:
            self.contig_ids = contig_ids
            self.seqrecords = self.get_seqrecords()
        elif seqrecords and not contig_ids:
            self.seqrecords = seqrecords
            self.contig_ids = self.get_contig_ids()
        else:
            self.seqrecords = seqrecords
            self.contig_ids = contig_ids
        self.nseqs = len(self.contig_ids)

    @property
    @lru_cache(maxsize=None)
    def assembly_seqs(self):
        """Retrieve all of the sequences from the provided `assembly`.

        Returns
        -------
        list
            list of str [seq, ...]

        """
        with open(self.assembly) as fh:
            return [seq for title, seq in SimpleFastaParser(fh)]

    def get_contig_ids(self):
        """Retrieve contig_ids from `self.seqrecords`.

            Returns
            -------
            set
                {contig_id, contig_id, ...}

            """
        return {seqrecord.id for seqrecord in self.seqrecords}

    def get_seqrecords(self):
        """Retrieve seqrecords from assembly contained in `self.contig_ids`.

        Returns
        -------
        list
            list of SeqIO [SeqRecords, ...]

        """
        return [
            seqrecord
            for seqrecord in SeqIO.parse(self.assembly, "fasta")
            if seqrecord.id in self.contig_ids
        ]

    @property
    @lru_cache(maxsize=None)
    def totalsize(self):
        """Get the total assembly size in bp.

        Returns
        -------
        int
            Sum of sequence lengths in `assembly`

        """
        return sum(len(rec) for rec in self.assembly_seqs)

    @property
    def size_pct(self):
        """Get the size of the `contigs` in reference to the `assembly` size.
        This is expressed as a percentage

        Returns
        -------
        float
            Size percentage of `contigs` covering the total `assembly` size.

        """
        return self.size / self.totalsize * 100

    @property
    def nallseqs(self):
        """Get the total number of sequences in `assembly`.

        Returns
        -------
        int
            Number of sequences contained in `assembly`

        """
        return len(self.assembly_seqs)

    @property
    def seqs_pct(self):
        """Get the percentage of the number of `contigs` compared to the number
        of sequences in the `assembly`.

        Returns
        -------
        float
            percentage of `contigs` corresponding to num. seqs. in `assembly`.

        """
        return self.nseqs / self.nallseqs * 100

    @property
    @lru_cache(maxsize=None)
    def size(self):
        """Get the summation of the lengths of `contigs` in bp.

        Returns
        -------
        int
            Summation of lengths of `contigs` in bp.

        """
        return sum(len(seq) for seq in self.seqrecords)

    def fragmentation_metric(self, quality_measure=0.50):
        """Describes the quality of assembled genomes that are fragmented in
        contigs of different length.

        For more information see:
            http://www.metagenomics.wiki/pdf/definition/assembly/n50

        Parameters
        ----------
        quality_measure : 0 < float < 1
            Description of parameter `quality_measure` (the default is .50).
            I.e. default measure is N50, but could use .1 for N10 or .9 for N90

        Returns
        -------
        int
            Minimum contig length to cover `quality_measure` of genome (i.e. percentile contig length)

        """
        target_size = self.size * quality_measure
        lengths = 0
        for length in sorted([len(record) for record in self.seqrecords], reverse=True):
            lengths += length
            if lengths > target_size:
                return length

    @property
    @lru_cache(maxsize=None)
    def length_weighted_gc(self):
        """Get the length weighted GC content of `seqrecords`.

        Returns
        -------
        float
            GC percentage weighted by contig length.

        """
        weights = [len(rec.seq) / self.size for rec in self.seqrecords]
        gc_content = [SeqUtils.GC(rec.seq) for rec in self.seqrecords]
        return np.average(a=gc_content, weights=weights)

    @property
    def nucl_orfs_exist(self):
        """Determine whether `nucl_orfs_fpath` exists.

        Returns
        -------
        bool
            True is returned if `nucl_orfs_fpath` is a valid path, otherwise False.

        """
        return self.prepared(self.nucl_orfs_fpath)

    @property
    def prot_orfs_exist(self):
        """Determine whether `prot_orfs_fpath` exists.

        Returns
        -------
        bool
            True is returned if `prot_orfs_fpath` is a valid path, otherwise False.

        """
        return self.prepared(self.prot_orfs_fpath)

    @property
    def nnucls(self):
        """Get the number of nucleotide ORFs in `nucl_orfs_fpath`.

        Returns
        -------
        int
            Number of nucleotide ORFs in `nucl_orfs_fpath`.

        """
        if not self.nucl_orfs_exist:
            return np.nan
        return len(self.get_orfs("nucl"))

    @property
    def nprots(self):
        """Get the number of amino-acid ORFs in `prot_orfs_fpath`.

        Returns
        -------
        int
            Number of amino-acid ORFs in `prot_orfs_fpath`.

        """
        if not self.prot_orfs_exist:
            return np.nan
        return len(self.get_orfs("prot"))

    def prepared(self, fpath):
        """Check whether `fpath` exists and is valid via checksum in checkpoints.

        Parameters
        ----------
        fpath : str
            </path/to/file>

        Returns
        -------
        bool
            Whether provided `fpath` has been generated and is valid

        """
        # COMBAK: Add checkpointing checksum check here
        if os.path.exists(fpath) and os.path.getsize(fpath) > 0:
            return True
        return False

    def get_orfs(self, orf_type="prot"):
        """Retrieve ORFs corresponding to MetaBin.

        Parameters
        ----------
        orf_type : str
            Type of ORF to retrieve (the default is 'prot'). Amino acid or nucleotide
            choices = ['prot','nucl']

        Returns
        -------
        list
            [SeqIO.SeqRecord, ...]

        Raises
        -------
        KeyError
            `orf_type` not in ['prot','nucl']
        FileNotFoundError
            Either `self.prot_orfs_fpath` or `self.nucl_orfs_fpath` not found.
            This will correspond to `orf_type`

        """
        orfs_fpaths = {"prot": self.prot_orfs_fpath, "nucl": self.nucl_orfs_fpath}
        if orf_type not in orfs_fpaths:
            raise KeyError(f"{orf_type} not in {orfs_fpaths}")
        orfs_fpath = orfs_fpaths[orf_type]
        if not self.prepared(orfs_fpath):
            raise FileNotFoundError(orfs_fpath)
        return prodigal.orf_records_from_contigs(
            contigs=self.contig_ids, fpath=self.prot_orfs_fpath
        )

    def write_orfs(self, fpath, orf_type="prot"):
        """Write `orf_type` ORFs from `contigs_ids` to `fpath`.

        Parameters
        ----------
        fpath : str
            </path/to/write/`orf_type`.orfs.fasta
        orf_type : str
            type of ORFs to write, choices=['prot','nucl']

        """
        records = self.get_orfs(orf_type=orf_type)
        SeqIO.write(records, fpath, "fasta")

    def describe(self, autometa_details=True):
        """Describe various statistics of the MetaBin instance.

        Parameters
        ----------
        autometa_details : bool
            If True, will output autometa related information (the default is True).

        """
        print(
            "Metabin Details\n"
            "________________________\n"
            f"Num. Contigs: {self.nseqs:,} / {self.nallseqs:,} ({self.seqs_pct:4.2f}% of total metagenome)\n"
            f"Num. Nucl. ORFs: {self.nnucls:,}\n"
            f"Num. Prot. ORFs: {self.nprots:,}\n"
            f"Size: {self.size:,}bp / {self.totalsize:,}bp ({self.size_pct:4.2f}% of total metagenome)\n"
            f"Length Weighted Avg. GC content: {self.length_weighted_gc:4.2f}%\n"
            "________________________\n"
        )
        if not autometa_details:
            return
        print(
            "Autometa Details\n"
            "________________________\n"
            f"Metagenome: {self.assembly}\n"
            f"Nucl. ORFs filepath: {self.nucl_orfs_fpath}\n"
            f"Prot. ORFs filepath: {self.prot_orfs_fpath}\n"
            f"Nucl. ORFs called: {self.nucl_orfs_exist}\n"
            f"Prot. ORFs called: {self.prot_orfs_exist}\n"
        )

    def subset_df(self, df):
        """Retrieve subset of provided `df` containing only `contig_ids`.

        Parameters
        ----------
        df : pd.DataFrame,pd.Series,str
            `df` may be a pandas Series or DataFrame, or a path to file.
            If a path is provided, `df` will be read as a tab-delimited table.

        Returns
        -------
        pd.DataFrame
            index='contig', cols=cols in provided `df`, subset by `contig_ids`.

        """
        if type(df) not in [pd.DataFrame, pd.Series]:
            raise TypeError(
                f"Unable to subset df. {type(df)} is not Series or DataFrame"
            )
        if type(df) is pd.DataFrame:
            return df[df.index.isin(self.contig_ids)]
        elif type(df) is str and self.prepared(df):
            df = pd.read_csv(df, sep="\t", index_col="contig")
        return df[df.index.isin(self.contig_ids)]


if __name__ == "__main__":
    print("MetaBin is part of Autometa and should not be run directly.")
    pass
