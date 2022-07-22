#!/usr/bin/env python

import logging

from abc import ABC, abstractmethod
from typing import Dict, Set, Tuple, List, Union, Iterable

import pandas as pd


logger = logging.getLogger(__name__)


class TaxonomyDatabase(ABC):

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

    def __init__(
        self,
        nodes: Dict[str, int],
    ) -> None:
        self.nodes = nodes

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
            Name of provided `taxid` if `taxid` is found in names.dmp else 'unclassified'

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
            Name of provided `taxid` if `taxid` is found in names.dmp else 'unclassified'

        """
        if not rank:
            return self.names.get(taxid, "unclassified")
        if rank not in set(TaxonomyDatabase.CANONICAL_RANKS):
            logger.warning(f"{rank} not in canonical ranks!")
            return "unclassified"
        ancestor_taxid = taxid
        while ancestor_taxid != 1:
            ancestor_rank = self.rank(ancestor_taxid)
            if ancestor_rank == rank:
                return self.names.get(ancestor_taxid, "unclassified")
            ancestor_taxid = self.parent(ancestor_taxid)
        # At this point we have not encountered a name for the taxid rank
        # so we will place this as unclassified.
        return "unclassified"

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
            rank name if taxid is found in nodes else "unclassified"

        """
        return self.nodes.get(taxid, {"rank": "unclassified"}).get("rank")

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
            df.fillna(value="unclassified", inplace=True)
        return df
