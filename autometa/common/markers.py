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

Autometa Marker class consisting of various methods to annotate sequences with
marker sets depending on sequence set taxonomy
"""

import logging
import os

import pandas as pd

from autometa.common.external import hmmer


BASE_DIR = os.path.dirname(os.path.dirname(__file__))
MARKERS_DIR = os.path.join(BASE_DIR,'databases','markers')

logger = logging.getLogger(__name__)


class Markers:
    """docstring for Autometa Markers class.

    Parameters
    ----------
    orfs_fpath : str
        Description of parameter `orfs_fpath`.
    kingdom : str, optional
        Description of parameter `kingdom` (the default is 'bacteria').
    dbdir : str, optional
        Description of parameter `dbdir` (the default is {MARKERS_DIR}).

    Attributes
    ----------
    hmmdb : str
        </path/to/hmmpressed/`dbdir`/`kingdom`.single_copy.hmm
    cutoffs : str
        </path/to/`dbdir`/`kingdom`.single_copy.cutoffs
    hmmscan_fn : str
        <`kingdom`.hmmscan.tsv>
    hmmscan_fp : str
        </path/to/`kingdom`.hmmscan.tsv>
    markers_fn : str
        <`kingdom`.markers.tsv>
    markers_fp : str
        </path/to/`kingdom`.markers.tsv>
    """
    def __init__(self, orfs_fpath, kingdom='bacteria', dbdir=MARKERS_DIR):
        self.orfs_fpath = os.path.realpath(orfs_fpath)
        self.kingdom = kingdom.lower()
        self.dbdir = dbdir
        self.hmmdb = os.path.join(self.dbdir, f'{self.kingdom}.single_copy.hmm')
        self.cutoffs = os.path.join(self.dbdir, f'{self.kingdom}.single_copy.cutoffs')
        self.hmmscan_fn = '.'.join([self.kingdom,'hmmscan.tsv'])
        self.hmmscan_fp = os.path.join(os.path.dirname(self.orfs_fpath),self.hmmscan_fn)
        self.markers_fn = '.'.join([self.kingdom,'markers.tsv'])
        self.markers_fp = os.path.join(os.path.dirname(self.orfs_fpath),self.markers_fn)

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
        NotImplementedError
            method not yet implemented

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
            </path/to/`kingdom`.markers.tsv>
        format : str
            * wide - index=contig, cols=[domain sacc,..] (default)
            * long - index=contig, cols=['sacc','count']
            * list - {contig:[sacc,...],...}
            * counts - {contig:len([sacc,...]), ...}

        Returns
        -------
        pd.DataFrame
            if 'wide' or 'long' else dict shape is (row x col)

        Raises
        -------
        FileNotFoundError
            Provided `fpath` does not exist
        ValueError
            Provided `format` is not in choices:
            choices = wide, long, list or counts

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
            # TODO: Write Marker specific AutometaException
            raise ValueError(err_msg)

    def get(self, format='wide', **kwargs):
        """Retrieve contigs' markers from markers database that pass cutoffs filter.

        Parameters
        ----------
        format : str, optional
            * wide - returns wide dataframe of contig PFAM counts (default)
            * long - returns long dataframe of contig PFAM counts
            * list - returns list of pfams for each contig
            * counts - returns count of pfams for each contig

        Returns
        -------
        pd.Dataframe or dict
            * wide - pd.DataFrame(index_col=contig, columns=[PFAM,...])
            * long - pd.DataFrame(index_col=contig, columns=['sacc','count'])
            * list - {contig:[pfam,pfam,...],contig:[...],...}
            * counts - {contig:count, contig:count,...}

        Raises
        -------
        ValueError
            Why the exception is raised.
        """
        if not self.searched:
            hmmer.hmmscan(
                orfs=self.orfs_fpath,
                hmmdb=self.hmmdb,
                outfpath=self.hmmscan_fp,
                **kwargs)
        if not self.found:
            hmmer.filter_markers(
                infpath=self.hmmscan_fp,
                outfpath=self.markers_fp,
                cutoffs=self.cutoffs,
                orfs=self.orfs_fpath)
        return Markers.load(fpath=self.markers_fp, format=format)

def main():
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser(description='Annotate ORFs with kingdom-marker information')
    parser.add_argument('orfs', help='</path/to/prot.orfs.faa>')
    parser.add_argument('kingdom', help='kingdom to search for markers',
        choices=['bacteria','archaea'], default='bacteria')
    parser.add_argument('--dbdir', help=f'</path/to/markers/dir> {MARKERS_DIR}',
        default=MARKERS_DIR)
    args = parser.parse_args()
    markers = Markers(orfs_fpath=args.orfs, kingdom=args.kingdom, dbdir=args.dbdir)
    markers.get()

if __name__ == '__main__':
    main()
