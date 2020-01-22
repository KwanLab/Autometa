#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Autometa Marker class consisting of various methods to annotate sequences with
marker sets depending on sequence set taxonomy
"""

import logging
import os

import pandas as pd

from autometa.common.external import hmmer


BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
MARKERS_DIR = os.path.join(BASE_DIR,'databases','markers')

logger = logging.getLogger(__name__)


class Markers:
    """docstring for Autometa Markers class.

    Parameters
    ----------
    orfs_fpath : type
        Description of parameter `orfs_fpath`.
    kingdom : type
        Description of parameter `kingdom` (the default is 'bacteria').
    dbdir : type
        Description of parameter `dbdir` (the default is MARKERS_DIR).

    Attributes
    ----------
    hmmdb : type
        Description of attribute `hmmdb`.
    cutoffs : type
        Description of attribute `cutoffs`.
    markers_fn : type
        Description of attribute `markers_fn`.
    markers_fp : type
        Description of attribute `markers_fp`.
    hmmscan_fn : type
        Description of attribute `hmmscan_fn`.
    hmmscan_fp : type
        Description of attribute `hmmscan_fp`.
    orfs_fpath
    kingdom
    dbdir

    """
    def __init__(self, orfs_fpath, kingdom='bacteria', dbdir=MARKERS_DIR):
        self.orfs_fpath = os.path.realpath(orfs_fpath)
        self.kingdom = kingdom.lower()
        self.dbdir = dbdir
        self.hmmdb = os.path.join(self.dbdir, f'{self.kingdom}.single_copy.hmm')
        self.cutoffs = os.path.join(self.dbdir, f'{self.kingdom}.single_copy.cutoffs')
        self.markers_fn = '.'.join([self.kingdom,'markers.tsv'])
        self.markers_fp = os.path.join(os.path.dirname(self.orfs_fpath),self.markers_fn)
        self.hmmscan_fn = '.'.join([self.kingdom,'hmmscan.tsv'])
        self.hmmscan_fp = os.path.join(os.path.dirname(self.orfs_fpath),self.hmmscan_fn)

    @property
    def searched(self):
        if os.path.exists(self.hmmscan_fp) and os.stat(self.hmmscan_fp).st_size > 0:
            return True
        return False

    @property
    def found(self):
        if os.path.exists(self.markers_fp) and os.stat(self.markers_fp).st_size > 0:
            return True
        return False

    @property
    def n_ucgs(self):
        return len(self.ucgs)

    def get_ucgs(self):
        """Retrieve Universally Conserved Markers Genes for sequences

        Parameters
        ----------
        kingdom : type
            Description of parameter `kingdom` (the default is self.kingdom).

        Returns
        -------
        type
            Description of returned object.

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        raise NotImplementedError

    def get_marker_set(self, ):
        raise NotImplementedError

    @staticmethod
    def load(fpath, format='wide'):
        """Read markers table into specified `format`.

        Parameters
        ----------
        fpath : str
            Description of parameter `fpath`.
        format : str
            Description of parameter `format` (the default is 'wide').

        Returns
        -------
        pd.DataFrame if 'wide' or 'long' else dict
        shape is (row x col)
            wide - index=contig, cols=[domain sacc,..]
            long - index=contig, cols=['sacc','count']
            list - {contig:[sacc,...],...}
            counts - {contig:len([sacc,...]), ...}

        Raises
        -------
        FileNotFoundError
            Provided `fpath` does not exist
        ValueError
            Provided `format` is not in choices:
            choices = ['wide','long','list','counts']

        """
        if not os.path.exists(fpath):
            raise FileNotFoundError(fpath)
        df = pd.read_csv(fpath, sep='\t', index_col='contig')
        grouped_df = df.groupby('contig')['sacc']
        if format == 'wide':
            return grouped_df.value_counts().unstack()
        elif format == 'long':
            return grouped_df.value_counts().reset_index(level=1, name='count')
        elif format == 'list':
            return {k:v.tolist() for k,v in list(grouped_df)}
        elif format == 'counts':
            return grouped_df.count().to_dict()
        else:
            params = ['wide','long','list','counts']
            err_msg = f'{format} is not a supported format.\n\tSupported formats: {params}'
            raise ValueError(err_msg)

    def get_markers(self, format='wide', **kwargs):
        """Retrieve contigs' markers from markers database that pass cutoffs filter.

        Parameters
        ----------
        format : str
            Description of parameter `format` (the default is 'wide').
            Choices include:
                - 'wide' returns wide dataframe of contig PFAM counts (default)
                - 'long' returns long dataframe of contig PFAM counts
                - 'list' returns list of pfams for each contig
                - 'counts' returns count of pfams for each contig

        Returns
        -------
        pd.Dataframe or dict
            wide - pd.DataFrame(index_col='contig', columns=[PFAM,...])
            long - pd.DataFrame(index_col='contig', columns=['sacc','count'])
            list - {contig:[pfam,pfam,...],contig:[...],...}
            counts - {contig:count, contig:count,...}

        Raises
        -------
        ValueError
            Why the exception is raised.
        """
        if not self.searched:
            hmmer.hmmscan(self.orfs_fpath, self.hmmdb, self.hmmscan_fp, **kwargs)
        if not self.found:
            hmmer.filter_markers(self.hmmscan_fp, self.markers_fp, self.cutoffs)
        return Markers.load(fpath=self.markers_fp, format=format)

def main(args):
    markers = Markers(orfs_fpath=args.orfs, kingdom=args.kingdom, dbdir=args.dbdir)
    markers.get_markers()

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('')
    parser.add_argument('orfs', help='</path/to/prot.orfs.faa>')
    parser.add_argument('kingdom', help='kingdom to search for markers',
        choices=['bacteria','archaea'], default='bacteria')
    parser.add_argument('--dbdir', help='</path/to/markers/dir>', default=MARKERS_DIR)
    args = parser.parse_args()
    main(args)
