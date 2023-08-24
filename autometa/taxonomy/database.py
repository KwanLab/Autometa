#!/usr/bin/env python

import logging
import string

from abc import ABC, abstractmethod
from typing import Dict, Set, Tuple, List, Union, Iterable
from autometa.common.exceptions import DatabaseOutOfSyncError

import pandas as pd


logger = logging.getLogger(__name__)


class TaxonomyDatabase(ABC):

    """
    TaxonomyDatabase Abstract Base Class

    Abstract methods

    1. parse_nodes(self)
    2. parse_names(self)
    3. parse_merged(self)
    4. parse_delnodes(self)
    5. convert_accessions_to_taxids(self)

    e.g.

    class GTDB(TaxonomyDatabase):
        def __init__(self, ...):
            self.nodes = self.parse_nodes()
            self.names = self.parse_names()
            self.merged = self.parse_merged()
            self.delnodes = self.parse_delnodes()
            ...
        def parse_nodes(self):
            ...
        def parse_nodes(self):
            ...
        def parse_merged(self):
            ...
        def parse_delnodes(self):
            ...
        def convert_accessions_to_taxids(self, accessions):
            ...

    Available methods (after aforementioned implementations):

    1. convert_taxid_dtype
    2. name
    3. rank
    4. parent
    5. lineage
    6. is_common_ancestor
    7. get_lineage_dataframe

    Available attributes:

    CANONICAL_RANKS
    UNCLASSIFIED
    """

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
    UNCLASSIFIED = "unclassified"

    @abstractmethod
    def parse_nodes(self) -> Dict[int, Dict[str, Union[str, int]]]:
        """
        Parse the `nodes.dmp` database and set to self.nodes.

        Returns
        -------
        dict
            {child_taxid:{'parent':parent_taxid,'rank':rank}, ...}
        """

    @abstractmethod
    def parse_names(self) -> Dict[int, str]:
        """
        Parses through the names.dmp in search of the given `taxid` and returns its name

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
            Name of provided `taxid` if `taxid` is found in names.dmp else TaxonomyDatabase.UNCLASSIFIED

        """

    @abstractmethod
    def parse_merged(self) -> Dict[int, int]:
        """
        Parses merged.dmp such that merged `taxid`s may be updated with their up-to-date `taxid`s

        Returns
        -------
        dict
            {old_taxid: new_taxid, ...}
        """

    @abstractmethod
    def parse_delnodes(self) -> Set[int]:
        """
        Parses delnodes.dmp such that deleted `taxid`s may be updated with their up-to-date `taxid`s

        Returns
        -------
        set
            {taxid, ...}
        """

    @abstractmethod
    def convert_accessions_to_taxids(
        self,
        accessions: Dict[str, Set[str]],
    ) -> Tuple[Dict[str, Set[int]], pd.DataFrame]:
        """
        Translates subject sequence ids to taxids

        Parameters
        ----------
        accessions : dict
            {qseqid: {sseqid, ...}, ...}

        Returns
        -------
        Tuple[Dict[str, Set[int]], pd.DataFrame]
            {qseqid: {taxid, taxid, ...}, ...}, index=range, cols=[qseqid, sseqid, raw_taxid, ..., cleaned_taxid]

        """

    def convert_taxid_dtype(self, taxid: int) -> int:
        """
        1. Converts the given `taxid` to an integer and checks whether it is positive.
        2. Checks whether `taxid` is present in both nodes.dmp and names.dmp.
        3a. If (2) is false, will check for corresponding `taxid` in merged.dmp and convert to this then redo (2).
        3b. If (2) is true, will return converted taxid.
        4. If (3a) is false will look for `taxid` in delnodes.dmp. If present will convert to root (taxid=1)

        Parameters
        ----------
        taxid : int
            identifier for a taxon in NCBI taxonomy databases - nodes.dmp, names.dmp or merged.dmp

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
        #  Step 1 Converts the given `taxid` to an integer and checks whether it is positive.
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
        #  Step 2: Check whether taxid is present in both nodes.dmp and names.dmp.
        if taxid not in self.names and taxid not in self.nodes:
            # Step 3a. Check for corresponding taxid in merged.dmp
            if taxid not in self.merged:
                # Step 4: look for taxid in delnodes.dmp. If present will convert to root (taxid=1)
                if taxid in self.delnodes:
                    # Assign deleted taxid to root...
                    if self.verbose:
                        logger.debug(
                            f"Found {taxid} in delnodes.dmp, converting to root (taxid=1)"
                        )
                    taxid = 1
                else:
                    err_message = f"Databases out of sync. {taxid} not in found in nodes.dmp, names.dmp, merged.dmp or delnodes.dmp"
                    raise DatabaseOutOfSyncError(err_message)
            else:
                # Step 3b. convert taxid from merged.
                if self.verbose:
                    logger.debug(
                        f"Converted {taxid} to {self.merged[taxid]} from merged.dmp"
                    )
                taxid = self.merged[taxid]
                if taxid not in self.names and taxid not in self.nodes:
                    # NOTE: Do not check delnodes.dmp here, as at this point it appears the databases are indeed out of sync.
                    # i.e. The taxid should either be in merged.dmp or in delnodes.dmp if not in nodes.dmp or names.dmp
                    err_message = f"Databases out of sync. Merged taxid ({taxid}) not found in nodes.dmp or names.dmp!"
                    raise DatabaseOutOfSyncError(err_message)
        return taxid

    def name(self, taxid: int, rank: str = None) -> str:
        """
        Parses through the names.dmp in search of the given `taxid` and returns its name.

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
            Name of provided `taxid` if `taxid` is found in names.dmp else TaxonomyDatabase.UNCLASSIFIED

        """
        try:
            taxid = self.convert_taxid_dtype(taxid)
        except DatabaseOutOfSyncError as err:
            logger.warning(err)
            taxid = 0
        if not rank:
            return self.names.get(taxid, TaxonomyDatabase.UNCLASSIFIED)
        if rank not in set(TaxonomyDatabase.CANONICAL_RANKS):
            logger.warning(f"{rank} not in canonical ranks!")
            return TaxonomyDatabase.UNCLASSIFIED
        ancestor_taxid = taxid
        while ancestor_taxid != 1:
            ancestor_rank = self.rank(ancestor_taxid)
            if ancestor_rank == rank:
                return self.names.get(ancestor_taxid, TaxonomyDatabase.UNCLASSIFIED)
            ancestor_taxid = self.parent(ancestor_taxid)
        # At this point we have not encountered a name for the taxid rank
        # so we will place this as unclassified.
        return TaxonomyDatabase.UNCLASSIFIED

    def rank(self, taxid: int) -> str:
        """
        Return the respective rank of provided `taxid`.

        Parameters
        ----------
        taxid : int
            `taxid` to retrieve rank from nodes

        Returns
        -------
        str
            rank name if taxid is found in nodes else autoattribute:: autometa.taxonomy.database.TaxonomyDatabase.UNCLASSIFIED

        """
        try:
            taxid = self.convert_taxid_dtype(taxid)
        except DatabaseOutOfSyncError as err:
            logger.warning(err)
            taxid = 0
        return self.nodes.get(taxid, {"rank": TaxonomyDatabase.UNCLASSIFIED}).get(
            "rank"
        )

    def parent(self, taxid: int) -> int:
        """
        Retrieve the parent taxid of provided `taxid`.

        Parameters
        ----------
        taxid : int
           child taxid to retrieve parent

        Returns
        -------
        int
            Parent taxid if found in nodes otherwise 1

        """
        try:
            taxid = self.convert_taxid_dtype(taxid)
        except DatabaseOutOfSyncError as err:
            logger.warning(err)
            taxid = 0
        return self.nodes.get(taxid, {"parent": 1}).get("parent")

    def lineage(
        self, taxid: int, canonical: bool = True
    ) -> List[Dict[str, Union[str, int]]]:
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
            if canonical and self.rank(taxid) not in TaxonomyDatabase.CANONICAL_RANKS:
                taxid = self.parent(taxid)
                continue
            lineage.append(
                {"taxid": taxid, "name": self.name(taxid), "rank": self.rank(taxid)}
            )
            taxid = self.parent(taxid)
        return lineage

    def is_common_ancestor(self, taxid_A: int, taxid_B: int) -> bool:
        """
        Determines whether the provided taxids have a non-root common ancestor

        Parameters
        ----------
        taxid_A : int
            taxid in taxonomy database
        taxid_B : int
            taxid in taxonomy database

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

    def get_lineage_dataframe(
        self, taxids: Iterable, fillna: bool = True
    ) -> pd.DataFrame:
        """
        Given an iterable of taxids generate a pandas DataFrame of their canonical
        lineages

        Parameters
        ----------
        taxids : iterable
            `taxids` whose lineage dataframe is being returned
        fillna : bool, optional
            Whether to fill the empty cells with TaxonomyDatabase.UNCLASSIFIED or not, default True

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
        canonical_ranks = [
            rank
            for rank in reversed(TaxonomyDatabase.CANONICAL_RANKS)
            if rank != "root"
        ]
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
            df.fillna(value=TaxonomyDatabase.UNCLASSIFIED, inplace=True)
        return df
