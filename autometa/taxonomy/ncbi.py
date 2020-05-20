#!/usr/bin/env python
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

File containing definition of the NCBI class and containing functions useful for handling NCBI taxonomy databases
"""


import logging
import os
import gzip
import subprocess
import sys
import pickle

import numpy as np
import pandas as pd

from tqdm import tqdm

from autometa.common.utilities import make_pickle, unpickle, timeit


logger = logging.getLogger(__name__)

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
NCBI_DIR = os.path.join(BASE_DIR, 'databases', 'ncbi')


class NCBI:
    """Taxonomy utilities for NCBI databases."""

    CANONICAL_RANKS = [
        'species',
        'genus',
        'family',
        'order',
        'class',
        'phylum',
        'superkingdom',
        'root'
    ]

    def __init__(self, dirpath, verbose=False):
        self.dirpath = dirpath
        self.verbose = verbose
        self.disable = not self.verbose
        self.names_fpath = os.path.join(self.dirpath, 'names.dmp')
        self.nodes_fpath = os.path.join(self.dirpath, 'nodes.dmp')
        self.merged_fpath = os.path.join(self.dirpath, 'merged.dmp')
        self.accession2taxid_fpath = os.path.join(
            self.dirpath, 'prot.accession2taxid')
        acc2taxid_gz = '.'.join([self.accession2taxid_fpath, 'gz'])
        if not os.path.exists(self.accession2taxid_fpath) and os.path.exists(acc2taxid_gz):
            self.accession2taxid_fpath = acc2taxid_gz
        self.nr_fpath = os.path.join(self.dirpath, 'nr.gz')
        nr_bname = os.path.splitext(os.path.basename(self.nr_fpath))[0]
        nr_dmnd_fname = '.'.join([nr_bname, 'dmnd'])
        nr_dmnd_fpath = os.path.join(self.dirpath, nr_dmnd_fname)
        if os.path.exists(nr_dmnd_fpath):
            self.nr_fpath = nr_dmnd_fpath
        if 'dmnd' not in os.path.basename(self.nr_fpath):
            # This check should probably be in the dependencies/databases file...
            logger.warning(
                f'DatabaseWarning: {self.nr_fpath} needs to be formatted for diamond!')
        self.nodes = self.parse_nodes()
        self.names = self.parse_names()
        self.merged = self.parse_merged()

    def __repr__(self):
        return str(self)

    def __str__(self):
        # Perhaps should place summary here of files that do or do not exist?
        return self.dirpath

    def name(self, taxid, rank=None):
        """
        Parses through the `names.dmp` in search of the given `taxid` and returns its name

        Parameters
        ----------
        taxid : int
            identifer for a taxon in the Taxonomy Database by NCBI
        rank : str, optional
            If  provided, will return `taxid` name at `rank`, by default None
            Must be a canonical rank, choices: species, genus, family, order, class, phylum, superkingdom
            Eg. self.name(562, 'genus') would return 'Escherichia', where 562 is the taxid for Escherichia coli

        Returns
        -------
        str or None
            Name of provided `taxid` if taxid is found in names.dmp else will return None

        """
        if not self.is_valid_taxid(taxid):
            raise ValueError(
                f'Taxid must be a positive integer! Given: {taxid}')
        if not rank:
            return self.names.get(taxid, 'unclassified')
        if rank not in set(NCBI.CANONICAL_RANKS):
            logger.warning(f'{rank} not in canonical ranks!')
            return None
        ancestor_taxid = taxid
        if ancestor_taxid == 1:
            logger.error(
                f'Taxid cannot be 1! Taxid 1 is the root node which links to itself')
            return ancestor_taxid
        while ancestor_taxid != 1:
            ancestor_rank = self.rank(ancestor_taxid)
            if ancestor_rank == rank:
                return self.names.get(ancestor_taxid, 'unclassified')
            ancestor_taxid = self.parent(ancestor_taxid)

    def lineage(self, taxid, canonical=True):
        """
        Returns the lineage of taxids encountered when traversing to root

        Parameters
        ----------
        taxid : int
            identifer for a taxon in the Taxonomy Database by NCBI
        canonical : bool, optional
            Returns the names of all the canonical ranks when True and an empty list when False, by default True

        Returns
        -------
        ordered list of dicts
            [{'taxid':taxid, 'rank':rank,'name':'name'}, ...]
        """
        lineage = []
        while taxid != 1:  # This prevents the 'root' to be added in the list
            if canonical and self.rank(taxid) not in NCBI.CANONICAL_RANKS:
                taxid = self.parent(taxid)
                continue
            lineage.append({
                'taxid': taxid,
                'name': self.name(taxid),
                'rank': self.rank(taxid)})
            taxid = self.parent(taxid)
        return lineage

    def get_lineage_dataframe(self, taxids, fillna=True):
        """
        Given an iterable of taxids generate a pandas DataFrame of their canonical
        lineages.

        Parameters
        ----------
        taxid : int
            identifer for a taxon in the Taxonomy Database by NCBI
        fillna : bool, optional
            Whether to fill the empty cells  with 'unclassified' or not, default True

        Returns
        -------
        pd.DataFrame
            index = taxid
            columns = [superkingdom,phylum,class,order,family,genus,species]

        NOTE
        ----

        If you would like to merge the returned DataFrame with another
        DataFrame. Let's say where you retrieved your taxids:

        Example
        -------

        .. code-block:: python

            merged_df = pd.merge(
                your_df,
                this_df,
                how='left',
                left_on=<taxid_column>,
                right_index=True)
        """
        canonical_ranks = [r for r in reversed(NCBI.CANONICAL_RANKS)]
        canonical_ranks.remove('root')
        taxids = list(set(taxids))
        ranked_taxids = {}
        for rank in canonical_ranks:
            for taxid in taxids:
                # if not self.is_valid_taxid(taxid):
                #     return False
                name = self.name(taxid, rank=rank)
                if taxid not in ranked_taxids:
                    ranked_taxids.update({taxid: {rank: name}})
                else:
                    ranked_taxids[taxid].update({rank: name})
        df = pd.DataFrame(ranked_taxids).transpose()
        df.index.name = 'taxid'
        if fillna:
            df.fillna(value='unclassified', inplace=True)
        return df

    def rank(self, taxid):
        """
        Takes a taxid and returns its rank

        Parameters
        ----------
        taxid : int
            identifer for a taxon in the Taxonomy Database by NCBI
            `taxid` to retrieve rank from `nodes`

        Returns
        -------
        bool or str
            Flase if taxid is not valid else rank name if taxid is found in
            nodes.dmp  else 'unclassified'
        """
        if not self.is_valid_taxid(taxid):
            raise ValueError(
                f'Taxid must be a positive integer! Given: {taxid}')
        return self.nodes.get(taxid, {'rank': 'unclassified'}).get('rank')

    def parent(self, taxid):
        """
        Retrieve the parent taxid of provided taxid

        Parameters
        ----------
        taxid : int
            identifer for a taxon in the Taxonomy Database by NCBI

        Returns
        -------
        bool or int
            Flase if taxid is not valid else parent taxid if taxid is found
            in nodes.dmp otehrwise 1
        """
        if not self.is_valid_taxid(taxid):
            raise ValueError(
                f'Taxid must be a positive integer! Given: {taxid}')
        return self.nodes.get(taxid, {'parent': 1}).get('parent')

    # @timeit

    def parse_names(self):
        """
        Parses through `names.dmp` database and loads taxids' which have scientific names

        Returns
        -------
        dict
            {taxid:name, ...}

        """
        if self.verbose:
            logger.debug(f'Processing names from {self.names_fpath}')
        names = {}
        fh = open(self.names_fpath)
        for line in tqdm(fh, disable=self.disable, desc='parsing names', leave=False):
            taxid, name, __, classification = line.strip(
                '\t|\n').split('\t|\t')[:4]
            taxid = int(taxid)
            name = name.lower()
            # Only add scientific name entries
            if classification == 'scientific name':
                names.update({taxid: name})
        fh.close()
        if self.verbose:
            logger.debug('names loaded')
        return names

    # @timeit
    def parse_nodes(self):
        """
        Parses through `nodes.dmp` database and returns a nested dictionary with child taxid as keys
        and corresponding parent taxid and rank as items

        Returns
        -------
        dict
            {child_taxid:{'parent':parent_taxid,'rank':rank}, ...}

        """
        if self.verbose:
            logger.debug(f'Processing nodes from {self.nodes_fpath}')
        fh = open(self.nodes_fpath)
        __ = fh.readline()  # root line
        nodes = {1: {'parent': 1, 'rank': 'root'}}
        for line in tqdm(fh, disable=self.disable, desc='parsing nodes', leave=False):
            child, parent, rank = line.split('\t|\t')[:3]
            parent, child = map(int, [parent, child])
            rank = rank.lower()
            nodes.update({child: {'parent': parent, 'rank': rank}})
        fh.close()
        if self.verbose:
            logger.debug('nodes loaded')
        return nodes

    def parse_merged(self):
        """
        Parses through `merged.dmp` databse and returns a dictionary having old and new taxid

        Returns
        -------
        dict
            {old_taxid: new_taxid, ...}

        """
        if self.verbose:
            logger.debug(f'Processing nodes from {self.merged_fpath}')
        fh = open(self.merged_fpath)
        merged = {}
        for line in tqdm(fh, disable=self.disable, desc='parsing merged', leave=False):
            old_taxid, new_taxid = [
                int(taxid) for taxid in line.strip('\t|\n').split('\t|\t')]
            merged.update({old_taxid: new_taxid})
        fh.close()
        if self.verbose:
            logger.debug('merged loaded')
        return merged

    def is_common_ancestor(self, taxid1, taxid2):
        """
        Determines whether the provided taxids have a non-root common ancestor

        Parameters
        ----------
        taxid1 : int
            identifer for a taxon in the Taxonomy Database by NCBI
        taxid2 : int
            identifer for a taxon in the Taxonomy Database by NCBI

        Returns
        -------
        boolean
            True if taxids share a common ancestor else False
        """
        taxid1_lineage = self.lineage(taxid1)
        taxid2_lineage = self.lineage(taxid2)
        taxid1_lineage_taxids = {taxids.get('taxid')
                                 for taxids in taxid1_lineage}
        taxid2_lineage_taxids = {taxids.get('taxid')
                                 for taxids in taxid2_lineage}
        common_ancestor = (taxid2_lineage_taxids.intersection(
            taxid1_lineage_taxids))
        common_ancestor.discard(1)  # This discards root if present
        return True if common_ancestor else False

    def is_valid_taxid(self, taxid):
        """
        Checks if the given taxid is a positive interger or not 

        Parameters
        ----------
        taxid : int
            identifer for a taxon in the Taxonomy Database by NCBI

        Returns
        -------
        boolean
            True if the 'taxid' is a positive interger  else False
        """

        # This check if an interger has been added as str, eg. "562"
        if isinstance(taxid, str):
            if taxid.isnumeric() or taxid.replace('.', '', 1).isnumeric():
                # checks if it is something like 12.0 vs. 12.9.
                if float(taxid).is_integer():
                    # See https://stackoverflow.com/a/47764450/12671809
                    taxid = int(float(taxid))
                else:
                    return False
            else:
                return False
        # `is_integer` returns error if only '12' is used, float(taxid) is needed
        if not float(taxid).is_integer() or taxid < 0:
            return False
        return True


def main():
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser(
        description='Autometa NCBI utilities class')
    parser.add_argument('ncbi', help='</path/to/ncbi/database/directory>')
    parser.add_argument(
        '--query-taxid',
        help='query taxid to retrieve its name, order, class, rank and complete lineage',
        default=1222,
        type=int
    )
    parser.add_argument('--verbose', help="add verbosity",
                        action='store_true', default=False)
    args = parser.parse_args()
    ncbi = NCBI(args.ncbi, args.verbose)
    query_taxid_parent = ncbi.parent(taxid=args.query_taxid)
    logger.info(
        f'{args.query_taxid} Name: {ncbi.name(taxid=args.query_taxid)}\n'
        f'{args.query_taxid} Order: {ncbi.name(taxid=args.query_taxid, rank="order")}\n'
        f'{args.query_taxid} Class: {ncbi.name(taxid=args.query_taxid, rank="class")}\n'
        f'{args.query_taxid} Rank: {ncbi.rank(taxid=args.query_taxid)}\n'
        f'{args.query_taxid} Lineage:\n{ncbi.lineage(taxid=args.query_taxid)}\n'
        f'{args.query_taxid} is_common_ancestor {query_taxid_parent}: '
        f'{ncbi.is_common_ancestor(taxid1=query_taxid_parent, taxid2=args.query_taxid)}\n'
    )


if __name__ == '__main__':
    main()
