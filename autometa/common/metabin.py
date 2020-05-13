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

from autometa.common.markers import Markers,MARKERS_DIR
from autometa.common import kmers
from autometa.common.utilities import timeit
from autometa.common.external import prodigal
from autometa.binning import recursive_dbscan


logger = logging.getLogger(__name__)


class MetaBin:
    """Autometa MetaBin class to manipulate/annotate contigs from a metagenome.

    Parameters
    ----------
    assembly : str
        </path/to/metagenome.fasta>
    contigs : list
        List of contigs to manipulate/annotate (must be contained in metagenome).
    outdir : str, optional
        </path/to/output/directory> (Default is the directory storing the `assembly`).

    Attributes
    ----------
    basename : str
        base name of `assembly`
    root : str
        root name of `assembly` (Will remove common extension like '.fasta')
    nucls_fname : str
        File name of `contigs` nucleotide ORFs
    prots_fname : str
        File name of `contigs` amino-acid ORFs
    nucl_orfs_fpath : str
        </path/to/nucleotide/`outdir`/`nucls_fname`>
    prot_orfs_fpath : str
        </path/to/amino/acid/`outdir`/`prots_fname`>
    nseqs : int
        Number of contigs in MetaBin

    """

    def __init__(self, assembly, contigs, outdir=None):
        self.assembly = os.path.realpath(assembly)
        self.basename = os.path.basename(self.assembly)
        self.root = os.path.splitext(self.basename)[0]
        self.outdir = os.path.realpath(outdir) if outdir else os.path.dirname(self.assembly)
        nucls_ext = 'orfs.fna'
        prots_ext = 'orfs.faa'
        self.nucls_fname = '.'.join([self.root, nucls_ext])
        self.prots_fname = '.'.join([self.root, prots_ext])
        self.nucl_orfs_fpath = os.path.join(self.outdir, self.nucls_fname)
        self.prot_orfs_fpath = os.path.join(self.outdir, self.prots_fname)
        self.contigs = contigs
        self.nseqs = len(self.contigs)

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
            return [seq for title,seq in SimpleFastaParser(fh)]

    @property
    @lru_cache(maxsize=None)
    def seqrecords(self):
        """Retrieve seqrecords from assembly contained in `self.contigs`.

        Returns
        -------
        list
            list of SeqIO [SeqRecords, ...]

        """
        return [seq for seq in SeqIO.parse(self.assembly, 'fasta')
            if seq.id in self.contigs]

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

    @property
    @lru_cache(maxsize=None)
    def length_weighted_gc(self):
        """Get the length weighted GC content of `contigs`.

        Returns
        -------
        float
            GC percentage weighted by contig length.

        """
        weights = [len(rec.seq)/self.size for rec in self.seqrecords]
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
        return len(self.get_orfs('nucl'))

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
        return len(self.get_orfs('prot'))

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

    def get_orfs(self, orf_type='prot'):
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
        orfs_fpaths = {'prot':self.prot_orfs_fpath, 'nucl':self.nucl_orfs_fpath}
        if orf_type not in orfs_fpaths:
            raise KeyError(f'{orf_type} not in {orfs_fpaths}')
        orfs_fpath = orfs_fpaths[orf_type]
        if not self.prepared(orfs_fpath):
            raise FileNotFoundError(orfs_fpath)
        return prodigal.orf_records_from_contigs(
            contigs=self.contigs,
            fpath=self.prot_orfs_fpath)

    def write_orfs(self, fpath, orf_type='prot'):
        """Write `orf_type` ORFs from `contigs` to `fpath`.

        Parameters
        ----------
        fpath : str
            </path/to/write/`orf_type`.orfs.fasta
        orf_type : str
            type of ORFs to write, choices=['prot','nucl']

        """
        records = self.get_orfs(orf_type=orf_type)
        SeqIO.write(records, fpath, 'fasta')

    def describe(self, autometa_details=True):
        """Describe various statistics of the MetaBin instance.

        Parameters
        ----------
        autometa_details : bool
            If True, will output autometa related information (the default is True).

        """
        print('Metabin Details\n'
            '________________________\n'
            f'Num. Contigs: {self.nseqs:,} / {self.nallseqs:,} ({self.seqs_pct:4.2f}% of total metagenome)\n'
            f'Num. Nucl. ORFs: {self.nnucls:,}\n'
            f'Num. Prot. ORFs: {self.nprots:,}\n'
            f'Size: {self.size:,}bp / {self.totalsize:,}bp ({self.size_pct:4.2f}% of total metagenome)\n'
            f'Length Weighted Avg. GC content: {self.length_weighted_gc:4.2f}%\n'
            '________________________\n')
        if not autometa_details:
            return
        print(
            'Autometa Details\n'
            '________________________\n'
            f'Metagenome: {self.assembly}\n'
            f'Nucl. ORFs filepath: {self.nucl_orfs_fpath}\n'
            f'Prot. ORFs filepath: {self.prot_orfs_fpath}\n'
            f'Nucl. ORFs called: {self.nucl_orfs_exist}\n'
            f'Prot. ORFs called: {self.prot_orfs_exist}\n')

    @timeit
    def get_binning(self, method='recursive_dbscan', **kwargs):
        """Retrieve binning results from provided `method`.

        Note: Most required arguments should be provided in `kwargs`, this is
        currently done to allow easy addition of other binning methods via `method` arg.

        Parameters
        ----------
        method : str
            Description of parameter `method` (the default is 'recursive_dbscan').
        **kwargs : dict
            Additional keyword arguments to be passed in respective to the provided
            `method`

            * kmers : str, </path/to/kmers.tsv>
            * embedded : str, </path/to/kmers.embedded.tsv>
            * do_pca : bool, Perform PCA prior to embedding
            * embedding_method : str, Embedding method to use choices=['sksne','bhsne','umap']
            * perplexity : float, perplexity setting for embedding method sksne/bhsne
            * coverage : str, </path/to/coverages.tsv>
            * taxonomy : str, </path/to/taxonomy.tsv>
            * domain : str, kindom to bin choices=['bacteria','archaea']
            * purity : float, purity cutoff to apply to bins
            * completeness : float, completeness cutoff to apply to bins
            * reverse : bool, If true will bin taxa from least to most specific rank

        Returns
        -------
        pd.DataFrame
            index=contig cols=['cluster','completeness','purity',...]

        Raises
        -------
        NotImplementedError
            Provided `method` is not yet implemented.

        """
        if method == 'recursive_dbscan':
            try:
                embedded_df = kmers.embed(
                    kmers=kwargs.get('kmers'),
                    embedded=kwargs.get('embedded'),
                    do_pca=kwargs.get('do_pca',True),
                    pca_dimensions=kwargs.get('pca_dims',50),
                    method=kwargs.get('embedding_method','bhsne'),
                    perplexity=kwargs.get('perplexity',30),
                )
            except ValueError as error:
                logger.exception(error)
            master_df = embedded_df
            coverage_fp = kwargs.get('coverage')
            if coverage_fp:
                coverage_df = pd.read_csv(coverage_fp, sep='\t', index_col='contig')
                master_df = pd.merge(
                    master_df,
                    coverage_df,
                    how='left',
                    left_index=True,
                    right_index=True)
            taxonomy_fp = kwargs.get('taxonomy')
            if taxonomy_fp:
                taxa_df = pd.read_csv(taxonomy_fp, sep='\t', index_col='contig')
                master_df = pd.merge(
                    left=master_df,
                    right=taxa_df,
                    how='left',
                    left_index=True,
                    right_index=True)
            master_df = self.subset_df(master_df)
            master_df = master_df.convert_dtypes()
            use_taxonomy = True if 'taxid' in master_df else False
            markers = self.markers(kwargs.get('domain','bacteria'))
            logger.info(f'Binning {kwargs.get("domain")} with {method}')
            return recursive_dbscan.binning(
                master=master_df,
                markers=markers,
                domain=kwargs.get('domain','bacteria'),
                completeness=kwargs.get('completeness',20.),
                purity=kwargs.get('purity',90.),
                taxonomy=use_taxonomy,
                method='DBSCAN',
                reverse=kwargs.get('reverse',True)
            )
        raise NotImplementedError(f'{method} not yet implemented')

    @timeit
    def markers(self, kingdom='bacteria', dbdir=MARKERS_DIR, force=False):
        f"""Retrieve Markers dataframe using orfs called from `orf_caller` and
        annotated belonging to provided `kingdom`.

        Parameters
        ----------
        kingdom : str, optional
            Domain specific markers to retrieve (the default is 'bacteria').
        dbdir : str, optional
            </path/to/markers/database/directory> (the default is {MARKERS_DIR}).
            Should contain pressed hmms and cutoffs table.
        force : bool, optional
            Will overwrite existing marker annotations (the default is {force}).

        Returns
        -------
        pd.DataFrame
            wide format - index_col='contig', columns=[PFAM,...]

        """

        logger.debug(f'Retrieving markers for {kingdom} kingdom')
        orfs_fp = os.path.join(self.outdir, f'{kingdom.lower()}.orfs.faa')
        if (not os.path.exists(orfs_fp)) or (os.path.exists(orfs_fp) and force):
            self.write_orfs(orfs_fp)
        markers = Markers(orfs_fp, kingdom=kingdom, dbdir=dbdir)
        return markers.get()

    def subset_df(self, df):
        """Retrieve subset of provided `df` containing only `contigs`.

        Parameters
        ----------
        df : pd.DataFrame,pd.Series,str
            `df` may be a pandas Series or DataFrame, or a path to file.
            If a path is provided, `df` will be read as a tab-delimited table.

        Returns
        -------
        pd.DataFrame
            index='contig', cols=cols in provided `df`, subset by `contigs`.

        """
        if type(df) not in [pd.DataFrame, pd.Series]:
            raise TypeError(f'Unable to subset df. {type(df)} is not Series or DataFrame')
        if type(df) is pd.DataFrame:
            return df[df.index.isin(self.contigs)]
        elif type(df) is str and self.prepared(df):
            df = pd.read_csv(df, sep='\t', index_col='contig')
        return df[df.index.isin(self.contigs)]

def main():
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser(description='Autometa MetaBin Class')
    parser.add_argument('--assembly', help='</path/to/metagenome.fasta>', required=True)
    parser.add_argument('--contigs', help='list of contigs in MetaBin',nargs='+', required=True)
    parser.add_argument('--domain', help='kingdom to use for binning', default='bacteria')
    parser.add_argument('--kmers', help='</path/to/kmers.tsv')
    parser.add_argument('--taxonomy', help='</path/to/taxonomy_vote.tsv')
    parser.add_argument('--coverage', help='</path/to/coverages.tsv')
    args = parser.parse_args()
    mag = MetaBin(args.assembly, args.contigs)

    mag.get_binning(
        method='recursive_dbscan',
        kmers=args.kmers,
        domain=args.domain,
        taxonomy=args.taxonomy,
        coverage=args.coverage)

if __name__ == '__main__':
    main()
