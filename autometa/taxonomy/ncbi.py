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
import string
import subprocess
import sys

import pandas as pd

from tqdm import tqdm

from autometa.common.utilities import make_pickle, unpickle, timeit
from autometa.common.exceptions import DatabaseOutOfSyncError


logger = logging.getLogger(__name__)

BASE_DIR = os.path.dirname(os.path.dirname(__file__))
NCBI_DIR = os.path.join(BASE_DIR, "databases", "ncbi")


class NCBI:
    """Taxonomy utilities for NCBI databases."""

    CANONICAL_RANKS = [
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "superkingdom",
        "root",
    ]

    def __init__(self, dirpath, verbose=False):
        """
        Instantiates the NCBI class

        Parameters
        ----------
        dirpath : str
            Path to the database directory
        verbose : bool, optional
            log progress to terminal, by default False
        """
        self.dirpath = dirpath
        self.verbose = verbose
        self.disable = not self.verbose
        self.names_fpath = os.path.join(self.dirpath, "names.dmp")
        self.nodes_fpath = os.path.join(self.dirpath, "nodes.dmp")
        self.merged_fpath = os.path.join(self.dirpath, "merged.dmp")
        self.accession2taxid_fpath = os.path.join(self.dirpath, "prot.accession2taxid")
        acc2taxid_gz = ".".join([self.accession2taxid_fpath, "gz"])
        if not os.path.exists(self.accession2taxid_fpath) and os.path.exists(
            acc2taxid_gz
        ):
            self.accession2taxid_fpath = acc2taxid_gz
        self.nr_fpath = os.path.join(self.dirpath, "nr.gz")
        nr_bname = os.path.splitext(os.path.basename(self.nr_fpath))[0]
        nr_dmnd_fname = ".".join([nr_bname, "dmnd"])
        nr_dmnd_fpath = os.path.join(self.dirpath, nr_dmnd_fname)
        if os.path.exists(nr_dmnd_fpath):
            self.nr_fpath = nr_dmnd_fpath
        if "dmnd" not in os.path.basename(self.nr_fpath):
            # This check should probably be in the dependencies/databases file...
            logger.warning(
                f"DatabaseWarning: {self.nr_fpath} needs to be formatted for diamond!"
            )
        self.nodes = self.parse_nodes()
        self.names = self.parse_names()
        self.merged = self.parse_merged()

    def __repr__(self):
        """
        Operator overloading to return the string representation of the class object

        Returns
        -------
        str
            String representation of the class object
        """
        return str(self)

    def __str__(self):
        """
        Operator overloading to return the directory path of the class object

        Returns
        -------
        str
            Directory path of the class object
        """
        # Perhaps should place summary here of files that do or do not exist?
        return self.dirpath

    def name(self, taxid, rank=None):
        """
        Parses through the names.dmp in search of the given `taxid` and returns its name. If the `taxid` is
        deprecated, suppressed, withdrawn from NCBI (basically old) the updated name will be retrieved

        Parameters
        ----------
        taxid : int
            `taxid` whose name is being returned
        rank : str, optional
            If  provided, will return `taxid` name at `rank`, by default None
            Must be a canonical rank, choices: species, genus, family, order, class, phylum, superkingdom
            Eg. self.name(562, 'genus') would return 'Escherichia', where 562 is the taxid for Escherichia coli

        Returns
        -------
        str
            Name of provided `taxid` if `taxid` is found in names.dmp else 'unclassified'

        Raises
        ------
        DatabaseOutOfSyncError
            NCBI databases nodes.dmp, names.dmp and merged.dmp are out of sync with each other
        """
        try:
            taxid = self.convert_taxid_dtype(taxid)
        except DatabaseOutOfSyncError as err:
            logger.warning(err)
            taxid = 0
        if not rank:
            return self.names.get(taxid, "unclassified")
        if rank not in set(NCBI.CANONICAL_RANKS):
            logger.warning(f"{rank} not in canonical ranks!")
            return "unclassified"
        ancestor_taxid = taxid
        while ancestor_taxid != 1:
            ancestor_rank = self.rank(ancestor_taxid)
            if ancestor_rank == rank:
                return self.names.get(ancestor_taxid, "unclassified")
            ancestor_taxid = self.parent(ancestor_taxid)
        # At this point, taxid must equal 1 (root)
        return self.names.get(ancestor_taxid, "unclassified")

    def lineage(self, taxid, canonical=True):
        """
        Returns the lineage of `taxids` encountered when traversing to root

        Parameters
        ----------
        taxid : int
            `taxid` in nodes.dmp, whose lineage is being returned
        canonical : bool, optional
            Lineage includes both canonical and non-canonical ranks when False, and only the canonical ranks when True
            Canonical ranks include : species, genus , family, order, class, phylum, superkingdom, root

        Returns
        -------
        ordered list of dicts
            [{'taxid':taxid, 'rank':rank,'name':name}, ...]
        """
        lineage = []
        while taxid != 1:
            if canonical and self.rank(taxid) not in NCBI.CANONICAL_RANKS:
                taxid = self.parent(taxid)
                continue
            lineage.append(
                {"taxid": taxid, "name": self.name(taxid), "rank": self.rank(taxid)}
            )
            taxid = self.parent(taxid)
        return lineage

    def get_lineage_dataframe(self, taxids, fillna=True):
        """
        Given an iterable of taxids generate a pandas DataFrame of their canonical
        lineages

        Parameters
        ----------
        taxids : iterable
            `taxids` whose lineage dataframe is being returned
        fillna : bool, optional
            Whether to fill the empty cells  with 'unclassified' or not, default True

        Returns
        -------
        pd.DataFrame
            index = taxid
            columns = [superkingdom,phylum,class,order,family,genus,species]

        Example
        -------

        If you would like to merge the returned DataFrame ('this_df') with another
        DataFrame ('your_df'). Let's say where you retrieved your taxids:

        .. code-block:: python

            merged_df = pd.merge(
                left=your_df,
                right=this_df,
                how='left',
                left_on=<taxid_column>,
                right_index=True)
        """
        canonical_ranks = [r for r in reversed(NCBI.CANONICAL_RANKS)]
        if "root" in canonical_ranks:
            canonical_ranks.remove("root")
        taxids = list(set(taxids))
        ranked_taxids = {}
        for rank in canonical_ranks:
            for taxid in taxids:
                name = self.name(taxid, rank=rank)
                if taxid not in ranked_taxids:
                    ranked_taxids.update({taxid: {rank: name}})
                else:
                    ranked_taxids[taxid].update({rank: name})
        df = pd.DataFrame(ranked_taxids).transpose()
        df.index.name = "taxid"
        if fillna:
            df.fillna(value="unclassified", inplace=True)
        return df

    def rank(self, taxid):
        """
        Return the respective rank of provided `taxid`. If the `taxid` is deprecated, suppressed,
        withdrawn from NCBI (basically old) the updated rank will be retrieved

        Parameters
        ----------
        taxid : int
            `taxid` to retrieve rank from nodes.dmp

        Returns
        -------
        str
            rank name if taxid is found in nodes.dmp else "unclassified"

        Raises
        ------
        DatabaseOutOfSyncError
            NCBI databases nodes.dmp, names.dmp and merged.dmp are out of sync with each other
        """
        try:
            taxid = self.convert_taxid_dtype(taxid)
        except DatabaseOutOfSyncError as err:
            logger.warning(err)
            taxid = 0
        return self.nodes.get(taxid, {"rank": "unclassified"}).get("rank")

    def parent(self, taxid):
        """
        Retrieve the parent taxid of provided `taxid`. If the `taxid` is deprecated, suppressed,
        withdrawn from NCBI (basically old) the updated parent will be retrieved

        Parameters
        ----------
        taxid : int
           child taxid to retrieve parent

        Returns
        -------
        int
            Parent taxid if found in nodes.dmp otherwise 1

        Raises
        ------
        DatabaseOutOfSyncError
            NCBI databases nodes.dmp, names.dmp and merged.dmp are out of sync with each other
        """
        try:
            taxid = self.convert_taxid_dtype(taxid)
        except DatabaseOutOfSyncError as err:
            logger.warning(err)
            taxid = 0
        return self.nodes.get(taxid, {"parent": 1}).get("parent")

    # @timeit
    def parse_names(self):
        """
        Parses through names.dmp database and loads taxids with scientific names

        Returns
        -------
        dict
            {taxid:name, ...}
        """
        if self.verbose:
            logger.debug(f"Processing names from {self.names_fpath}")
        names = {}
        fh = open(self.names_fpath)
        for line in tqdm(fh, disable=self.disable, desc="parsing names", leave=False):
            taxid, name, __, classification = line.strip("\t|\n").split("\t|\t")[:4]
            taxid = int(taxid)
            name = name.lower()
            # Only add scientific name entries
            if classification == "scientific name":
                names.update({taxid: name})
        fh.close()
        if self.verbose:
            logger.debug("names loaded")
        return names

    # @timeit
    def parse_nodes(self):
        """
        Parse the `nodes.dmp` database to be used later by :func:`autometa.taxonomy.ncbi.NCBI.parent`, :func:`autometa.taxonomy.ncbi.NCBI.rank`
        Note: This is performed when a new NCBI class instance is constructed

        Returns
        -------
        dict
            {child_taxid:{'parent':parent_taxid,'rank':rank}, ...}
        """
        if self.verbose:
            logger.debug(f"Processing nodes from {self.nodes_fpath}")
        fh = open(self.nodes_fpath)
        __ = fh.readline()  # root line
        nodes = {1: {"parent": 1, "rank": "root"}}
        for line in tqdm(fh, disable=self.disable, desc="parsing nodes", leave=False):
            child, parent, rank = line.split("\t|\t")[:3]
            parent, child = map(int, [parent, child])
            rank = rank.lower()
            nodes.update({child: {"parent": parent, "rank": rank}})
        fh.close()
        if self.verbose:
            logger.debug("nodes loaded")
        return nodes

    def parse_merged(self):
        """
        Parse the merged.dmp database
        Note: This is performed when a new NCBI class instance is constructed

        Returns
        -------
        dict
            {old_taxid: new_taxid, ...}
        """
        if self.verbose:
            logger.debug(f"Processing nodes from {self.merged_fpath}")
        fh = open(self.merged_fpath)
        merged = {}
        for line in tqdm(fh, disable=self.disable, desc="parsing merged", leave=False):
            old_taxid, new_taxid = [
                int(taxid) for taxid in line.strip("\t|\n").split("\t|\t")
            ]
            merged.update({old_taxid: new_taxid})
        fh.close()
        if self.verbose:
            logger.debug("merged loaded")
        return merged

    def is_common_ancestor(self, taxid_A, taxid_B):
        """
        Determines whether the provided taxids have a non-root common ancestor

        Parameters
        ----------
        taxid_A : int
            taxid in NCBI taxonomy databases - nodes.dmp, names.dmp or merged.dmp
        taxid_B : int
            taxid in NCBI taxonomy databases - nodes.dmp, names.dmp or merged.dmp

        Returns
        -------
        boolean
            True if taxids share a common ancestor else False
        """
        lineage_a_taxids = {ancestor.get("taxid") for ancestor in self.lineage(taxid_A)}
        lineage_b_taxids = {ancestor.get("taxid") for ancestor in self.lineage(taxid_B)}
        common_ancestor = lineage_b_taxids.intersection(lineage_a_taxids)
        common_ancestor.discard(1)  # This discards root
        return True if common_ancestor else False

    def convert_taxid_dtype(self, taxid):
        """
        1. Converts the given `taxid` to an integer and checks whether it is positive.
        2. Checks whether `taxid` is present in both nodes.dmp and names.dmp.
        3a. If (2) is false, will check for corresponding `taxid` in merged.dmp and convert to this then redo (2).
        3b. If (2) is true, will return converted taxid.

        Parameters
        ----------
        taxid : int
            identifer for a taxon in NCBI taxonomy databases - nodes.dmp, names.dmp or merged.dmp

        Returns
        -------
        int
            `taxid` if the `taxid` is a positive integer and present in either nodes.dmp or names.dmp or
            taxid recovered from merged.dmp

        Raises
        ------
        ValueError
            Provided `taxid` is not a positive integer
        DatabaseOutOfSyncError
            NCBI databases nodes.dmp, names.dmp and merged.dmp are out of sync with each other
        """
        # Checking taxid instance format
        # This checks if an integer has been added as str, eg. "562"
        if isinstance(taxid, str):
            invalid_chars = set(string.punctuation + string.ascii_letters)
            invalid_chars.discard(".")
            if set(taxid).intersection(invalid_chars) or taxid.count(".") > 1:
                raise ValueError(f"taxid contains invalid character(s)! Given: {taxid}")
            taxid = float(taxid)
        # a boolean check is needed as they will evaluate silently to 0 or 1 when cast as ints. FALSE=0, TRUE=1
        # float(taxid).is_integer() checks if it is something like 12.0 vs. 12.9
        # is_integer only takes float as input else raises error, thus isinstance( ,float) is used before it to make sure a float is being passed
        if isinstance(taxid, bool) or (
            isinstance(taxid, float) and not taxid.is_integer()
        ):
            raise ValueError(f"taxid must be an integer! Given: {taxid}")
        taxid = int(taxid)
        if taxid <= 0:
            raise ValueError(f"Taxid must be a positive integer! Given: {taxid}")
        # Checking databases
        if taxid not in self.names and taxid not in self.nodes:
            if taxid not in self.merged:
                err_message = f"Databases out of sync. {taxid} not in found in nodes.dmp, names.dmp or merged.dmp."
                raise DatabaseOutOfSyncError(err_message)
            else:
                taxid = self.merged[taxid]
                if taxid not in self.names and taxid not in self.nodes:
                    err_message = f"Databases out of sync. Merged taxid ({taxid}) not found in nodes.dmp or names.dmp!"
                    raise DatabaseOutOfSyncError(err_message)
        return taxid


if __name__ == "__main__":
    print("file containing Autometa NCBI utilities class")
    import sys

    sys.exit(1)
