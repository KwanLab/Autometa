#!/usr/bin/env python
"""
Utilities file containing functions useful for handling NCBI taxonomy databases
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

from autometa.common.utilities import make_pickle, unpickle,timeit


logger = logging.getLogger(__name__)

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
NCBI_DIR = os.path.join(BASE_DIR,'databases','ncbi')

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
        self.accession2taxid_fpath = os.path.join(self.dirpath, 'prot.accession2taxid')
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
            logger.warning(f'DatabaseWarning: {self.nr_fpath} needs to be formatted for diamond!')
        self.nodes = self.parse_nodes()
        self.names = self.parse_names()
        self.merged = self.parse_merged()

    def __repr__(self):
        return str(self)

    def __str__(self):
        # Perhaps should place summary here of files that do or do not exist?
        return self.dirpath

    def name(self, taxid, rank=None):
        """return `taxid` name

        Parameters
        ----------
        taxid : int
            Description of parameter `taxid`.
        rank : str
            If  provided, will return `taxid` name at `rank` (the default is None).
            Must be a canonical rank.
            Choices: species, genus, family, order, class, phylum, superkingdom

        Returns
        -------
        str
            Name of provided `taxid`.

        Raises
        -------
        ValueError
            Provided `taxid` must be an integer

        """
        # If taxid is not found in names.dmp, will return None
        try:
            taxid = int(taxid)
        except ValueError as err:
            logger.error(f'Taxid must be an integer! {taxid} type --> {type(taxid)}')
            return None
        if not rank:
            return self.names.get(taxid, 'unclassified')
        if rank not in set(NCBI.CANONICAL_RANKS):
            logger.warning(f'{rank} not in canonical ranks!')
            return None
        ancestor_taxid = taxid
        while ancestor_taxid != 1:
            ancestor_rank = self.rank(ancestor_taxid)
            if ancestor_rank == rank:
                return self.names.get(ancestor_taxid, 'unclassified')
            ancestor_taxid = self.parent(ancestor_taxid)

    def lineage(self, taxid, canonical=True):
        """Returns the lineage of taxids encountered when traversing to root

        Parameters
        ----------
        taxid : int
            `taxid` in nodes.dmp

        Returns
        -------
        ordered list of dicts
            [{'taxid':taxid, rank':rank,'name':'name'}, ...]

        Raises
        -------
        ValueError
            Provided `taxid` is not an integer
        """
        try:
            taxid = int(taxid)
        except ValueError as err:
            logger.error(f'Taxid must be an integer! {taxid} type --> {type(taxid)}')
            return None
        lineage = []
        while taxid != 1:
            if canonical and self.rank(taxid) not in NCBI.CANONICAL_RANKS:
                taxid = self.parent(taxid)
                continue
            lineage.append({
                'taxid':taxid,
                'name':self.name(taxid),
                'rank':self.rank(taxid)})
            taxid = self.parent(taxid)
        return lineage

    def get_lineage_dataframe(self, taxids, fillna=True):
        """Given an iterable of taxids generate a pandas DataFrame of their canonical
        lineages.

        Returns
        -------
        pd.DataFrame
            index = taxid
            columns = [superkingdom,phylum,class,order,family,genus,species]

        NOTE: If you would like to merge the returned DataFrame with another
        DataFrame... Let's say where you retrieved your taxids...

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
        taxids_ = {}
        for rank in canonical_ranks:
            for taxid in taxids:
                name = self.name(taxid, rank=rank)
                if taxid not in taxids_:
                    taxids_.update({taxid:{rank:name}})
                else:
                    taxids_[taxid].update({rank:name})
        df = pd.DataFrame(taxids_).transpose()
        df.index.name = 'taxid'
        if fillna:
            df.fillna(value='unclassified', inplace=True)
        return df

    def rank(self, taxid):
        """Short summary.

        Parameters
        ----------
        taxid : int
            `taxid` to retrieve rank from `nodes`.

        Returns
        -------
        str
            rank name

        Raises
        -------
        ValueError
            Provided `taxid` is not an integer

        """
        # If taxid is not found in nodes.dmp, will return 'unclassified'
        try:
            taxid = int(taxid)
        except ValueError as err:
            logger.error(f'Taxid must be an integer! {taxid} type --> {type(taxid)}')
            return None
        return self.nodes.get(taxid,{'rank':'unclassified'}).get('rank')

    def parent(self, taxid):
        """Retrieve the parent taxid of provided taxid

        Parameters
        ----------
        taxid : int
            Description of parameter `taxid`.

        Returns
        -------
        int
            parent taxid
        """
        # If taxid is not found in nodes.dmp, will return 1
        try:
            taxid = int(taxid)
        except ValueError as err:
            logger.error(f'Taxid must be an integer! {taxid} type --> {type(taxid)}')
            return None
        return self.nodes.get(taxid,{'parent':1}).get('parent')

    # @timeit
    def parse_names(self):
        """Short summary.

        Parameters
        ----------
        names_fpath : str
            </path/to/names.dmp> `names_fpath`.
        verbose : boolean
            prints progress to terminal (the default is False).

        Returns
        -------
        dict
            {taxid:name, ...}

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        if self.verbose:
            logger.debug(f'Processing names from {self.names_fpath}')
        names = {}
        fh = open(self.names_fpath)
        for line in tqdm(fh, disable=self.disable, desc='parsing names', leave=False):
            taxid, name, __, classification = line.strip('\t|\n').split('\t|\t')[:4]
            taxid = int(taxid)
            name = name.lower()
            # Only add scientific name entries
            is_scientific = classification == 'scientific name'
            if is_scientific:
                names.update({taxid:name})
        fh.close()
        if self.verbose:
            logger.debug('names loaded')
        return names

    # @timeit
    def parse_nodes(self):
        """Short summary.

        Parameters
        ----------
        verbose : boolean
            prints progress to terminal (the default is False).

        Returns
        -------
        dict
            {child_taxid:{'parent':parent_taxid,'rank':rank}, ...}

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        if self.verbose:
            logger.debug(f'Processing nodes from {self.nodes_fpath}')
        fh = open(self.nodes_fpath)
        __ = fh.readline() # root line
        nodes = {1:{'parent':1, 'rank':'root'}}
        for line in tqdm(fh, disable=self.disable, desc='parsing nodes', leave=False):
            child, parent, rank = line.split('\t|\t')[:3]
            parent, child = map(int,[parent, child])
            rank = rank.lower()
            nodes.update({child:{'parent':parent,'rank':rank}})
        fh.close()
        if self.verbose:
            logger.debug('nodes loaded')
        return nodes

    def parse_merged(self):
        """Short summary.

        Parameters
        ----------
        verbose : boolean
            prints progress to terminal (the default is False).

        Returns
        -------
        dict
            {child_taxid:{'parent':parent_taxid,'rank':rank}, ...}

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        if self.verbose:
            logger.debug(f'Processing nodes from {self.merged_fpath}')
        fh = open(self.merged_fpath)
        merged = {}
        for line in tqdm(fh, disable=self.disable, desc='parsing merged', leave=False):
            old_taxid, new_taxid = [int(taxid) for taxid in line.strip('\t|\n').split('\t|\t')]
            merged.update({old_taxid:new_taxid})
        fh.close()
        if self.verbose:
            logger.debug('merged loaded')
        return merged

    def is_common_ancestor(self, parent_taxid, child_taxid):
        """Determines whether the provided taxids have a non-root common ancestor

        Parameters
        ----------
        parent_taxid : int
            Description of parameter `parent_taxid`.
        child_taxid : int
            Description of parameter `child_taxid`.

        Returns
        -------
        boolean
            True if taxids share a common ancestor and False otherwise
        """
        ancestor_taxid = child_taxid
        while ancestor_taxid != 1:
            if parent_taxid == ancestor_taxid:
                return True
            ancestor_taxid = self.parent(ancestor_taxid)
        return False

def main(args):
    ncbi = NCBI(args.ncbi, args.verbose)
    query_taxid_parent = ncbi.parent(taxid=args.query_taxid)
    logger.info(
        f'{args.query_taxid} Name: {ncbi.name(taxid=args.query_taxid)}\n'
        f'{args.query_taxid} Order: {ncbi.name(taxid=args.query_taxid, rank="order")}\n'
        f'{args.query_taxid} Class: {ncbi.name(taxid=args.query_taxid, rank="class")}\n'
        f'{args.query_taxid} Rank: {ncbi.rank(taxid=args.query_taxid)}\n'
        f'{args.query_taxid} Lineage:\n{ncbi.lineage(taxid=args.query_taxid)}\n'
        f'{args.query_taxid} is_common_ancestor {query_taxid_parent}: '
        f'{ncbi.is_common_ancestor(parent_taxid=query_taxid_parent, child_taxid=args.query_taxid)}\n'
    )

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Autometa NCBI utilities class')
    parser.add_argument('ncbi', help='</path/to/ncbi/database/directory>')
    parser.add_argument(
        '--query-taxid',
        help='query taxid to test NCBI class functionality',
        default=1222,
        type=int,
    )
    parser.add_argument('--verbose', help="add verbosity", action='store_true', default=True)
    args = parser.parse_args()
    main(args)
