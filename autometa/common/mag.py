#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
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

Autometa Bin Class
"""


import logging
import os

import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio import SeqUtils

from autometa.common.markers import Markers,MARKERS_DIR
from autometa.common import kmers
from autometa.common.utilities import timeit
from autometa.common.external import prodigal
from autometa.binning import recursive_dbscan


logger = logging.getLogger(__name__)


class MAG:
    """docstring for Autometa MAG class."""

    def __init__(self, assembly, contigs, outdir=None):
        self.assembly = os.path.realpath(assembly)
        self.outdir = os.path.realpath(outdir) if outdir else os.path.dirname(self.assembly)
        self.basename = os.path.basename(self.assembly)
        self.assembly_name = self.basename.split('.')[0]
        nucls_ext = 'orfs.fna'
        prots_ext = 'orfs.faa'
        self.root = os.path.splitext(self.basename)[0]
        self.nucls_fname = '.'.join([self.root, nucls_ext])
        self.prots_fname = '.'.join([self.root, prots_ext])
        self.nucl_orfs_fpath = os.path.join(self.outdir, self.nucls_fname)
        self.prot_orfs_fpath = os.path.join(self.outdir, self.prots_fname)
        self.contigs = contigs
        self.nseqs = len(self.contigs)

    @property
    def totalsize(self):
        return sum([len(rec) for rec in self.get_seqs(all=True)])

    @property
    def size_pct(self):
        return self.size / self.totalsize * 100

    @property
    def nallseqs(self):
        return len(self.get_seqs(all=True))

    @property
    def seqs_pct(self):
        return self.nseqs / self.nallseqs * 100

    @property
    def size(self):
        return sum(len(seq) for seq in self.get_seqs())

    @property
    def gc_content(self):
        return np.mean([SeqUtils.GC(rec.seq) for rec in self.get_seqs()])

    @property
    def nucl_orfs_exist(self):
        return self.prepared(self.nucl_orfs_fpath)

    @property
    def prot_orfs_exist(self):
        return self.prepared(self.prot_orfs_fpath)

    @property
    def nnucls(self):
        if not self.nucl_orfs_exist:
            return np.nan
        return len(self.get_orfs('nucl'))

    @property
    def nprots(self):
        if not self.prot_orfs_exist:
            return np.nan
        return len(self.get_orfs('prot'))

    def prepared(self, fpath):
        if os.path.exists(fpath) and os.stat(fpath).st_size > 0:
            return True
        return False

    def get_orfs(self, orf_type='prot'):
        """Retrieve ORFs corresponding to MAG.

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
        FileNotFoundError
            Why the exception is raised.

        """
        orfs_fpaths = {'prot':self.prot_orfs_fpath, 'nucl':self.nucl_orfs_fpath}
        orfs_fpath = orfs_fpaths.get(orf_type, 'prot')
        if not self.prepared(orfs_fpath):
            raise FileNotFoundError(orfs_fpath)
        return prodigal.orf_records_from_contigs(
            contigs=self.contigs,
            fpath=self.prot_orfs_fpath)

    def write_orfs(self, fpath, orf_type='prot'):
        records = self.get_orfs(orf_type=orf_type)
        SeqIO.write(records, fpath, 'fasta')

    def get_seqs(self, all=False):
        """Retrieve sequences from assembly.

        Parameters
        ----------
        all : bool
            Gets all sequences from assembly if True else sequences for MAG
            (the default is False).

        Returns
        -------
        list
            list of SeqIO [SeqRecords, ...]

        """
        if all:
            return [seq for seq in SeqIO.parse(self.assembly, 'fasta')]
        return [seq for seq in SeqIO.parse(self.assembly, 'fasta')
            if seq.id in self.contigs]

    def describe(self, autometa_details=True):
        print(f'''
M.A.G. Details
________________________
Num. Contigs: {self.nseqs:,} / {self.nallseqs:,} ({self.seqs_pct:4.2f}% of total metagenome)
Num. Nucl. ORFs: {self.nnucls:,}
Num. Prot. ORFs: {self.nprots:,}
Size: {self.size:,}bp / {self.totalsize:,}bp ({self.size_pct:4.2f}% of total metagenome)
Mean GC content: {self.gc_content:4.2f}%
________________________
''')
        if autometa_details:
            print(f'''
            Autometa Details
            ________________________
            Metagenome: {self.assembly}
            Nucl. ORFs filepath: {self.nucl_orfs_fpath}
            Prot. ORFs filepath: {self.prot_orfs_fpath}
            Nucl. ORFs called: {self.nucl_orfs_exist}
            Prot. ORFs called: {self.prot_orfs_exist}''')

    def split_nucleotides(self, ):
        raise NotImplementedError()

    @timeit
    def get_binning(self, method='recursive_dbscan', **kwargs):
        """Retrieve binning results from provided `method`.

        Parameters
        ----------
        method : str
            Description of parameter `method` (the default is 'recursive_dbscan').
        **kwargs : dict
            Additional keyword arguments to be passed in respective to the provided
            `method`

        Returns
        -------
        pd.DataFrame
            index=contig cols=['cluster','completeness','purity',...]

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        if method == 'recursive_dbscan':
            try:
                embedded_df = kmers.embed(
                    kmers=kwargs.get('kmers'),
                    embedded=kwargs.get('embedded'),
                    do_pca=kwargs.get('do_pca',True),
                    pca_dimensions=kwargs.get('pca_dims',50),
                    method=kwargs.get('embedding_method','UMAP'),
                    perplexity=kwargs.get('perplexity',30),
                )
            except ValueError as error:
                logger.exception(error)
            master_df = embedded_df
            coverage_df = kwargs.get('coverage')
            if coverage_df is not None or not coverage_df.empty:
                master_df = pd.merge(
                    master_df,
                    coverage_df,
                    how='left',
                    left_index=True,
                    right_index=True)
            if kwargs.get('taxonomy'):
                taxa_df = pd.read_csv(kwargs.get('taxonomy'), sep='\t', index_col='contig')
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

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """

        logger.debug(f'Retrieving markers for {kingdom} kingdom')
        orfs_fp = os.path.join(self.outdir, f'{kingdom.lower()}.orfs.faa')
        if (not os.path.exists(orfs_fp)) or (os.path.exists(orfs_fp) and force):
            self.write_orfs(orfs_fp)
        markers = Markers(orfs_fp, kingdom=kingdom, dbdir=dbdir)
        return markers.get()

    def subset_df(self, df):
        if type(df) not in [pd.DataFrame, pd.Series]:
            raise TypeError(f'Unable to subset df. {type(df)} is not Series or DataFrame')
        if type(df) is pd.DataFrame:
            return df[df.index.isin(self.contigs)]
        elif type(df) is str and self.prepared(df):
            df = pd.read_csv(df, sep='\t', index_col='contig')
        return df[df.index.isin(self.contigs)]

def main(args):
    mag = MAG(args.assembly, args.contigs)

    mag.get_binning(
        method='recursive_dbscan',
        kmers=args.kmers,
        domain=args.domain,
        taxonomy=args.taxonomy,
        coverage=args.coverage)
    # mag.split_nucleotides()

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Autometa Bin Class')
    parser.add_argument('--assembly', help='</path/to/metagenome.fasta>', required=True)
    parser.add_argument('--contigs', help='list of contigs in MAG',nargs='+', required=True)
    parser.add_argument('--domain', help='kingdom to use for binning', default='bacteria')
    parser.add_argument('--kmers', help='</path/to/kmers.tsv')
    parser.add_argument('--taxonomy', help='</path/to/taxonomy_vote.tsv')
    parser.add_argument('--coverage', help='</path/to/coverages.tsv')
    args = parser.parse_args()
    main(args)
