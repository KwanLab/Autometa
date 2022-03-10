#!/usr/bin/env python
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

This script contains the LCA class containing methods to determine the
Lowest Common Ancestor given a tab-delimited BLAST table, fasta file, or
iterable of SeqRecords.

Note: LCA will assume the BLAST results table is in output format 6.
"""


import functools
import logging
import os

from typing import Dict, Set, Tuple
from itertools import chain

import pandas as pd
import numpy as np

from tqdm import tqdm
from pickle import UnpicklingError

from autometa.config.environ import get_versions
from autometa.taxonomy.ncbi import NCBI, NCBI_DIR
from autometa.common.utilities import make_pickle, unpickle, file_length
from autometa.common.external import diamond
from autometa.common.external import prodigal


logger = logging.getLogger(__name__)


class LCA(NCBI):
    """LCA class containing methods to retrieve the Lowest Common Ancestor.

    LCAs may be computed given taxids, a fasta or BLAST results.

    Parameters
    ----------
    dbdir : str
        Path to directory containing files: nodes.dmp, names.dmp, merged.dmp, prot.accession2taxid.gz
    outdir : str
        Output directory path to to serialize intermediate files to disk for later lookup
    verbose : bool, optional
        Add verbosity to logging stream (the default is False).

    Attributes
    ----------
    disable : bool
        Opposite of verbose. Used to disable `tqdm` module.
    tour_fp : str
        </path/to/serialized/file/eulerian/tour.pkl.gz>
    tour : list
        Eulerian tour containing branches and leaves information from tree traversal.
    level_fp : str
        </path/to/serialized/file/level.pkl.gz>
    level : list
        Lengths from root corresponding to `tour` during tree traversal.
    occurrence_fp : str
        </path/to/serialized/file/level.pkl.gz>
    occurrence : dict
        Contains first occurrence of each taxid while traversing tree (index in `tour`).
        e.g. {taxid:index, taxid: index, ...}
    sparse_fp : str
        </path/to/serialized/file/sparse.pkl.gz>
    sparse : numpy.ndarray
        Precomputed LCA values corresponding to `tour`,`level` and `occurrence`.
    lca_prepared : bool
        Whether LCA internals have been computed (e.g. `tour`,`level`,`occurrence`,`sparse`).

    """

    def __init__(self, dbdir: str, verbose: bool = False, cache: str = None):
        super().__init__(dbdir, verbose=verbose)
        self.verbose = verbose
        self.dbdir = dbdir
        self.cache = cache
        self.disable = False if verbose else True
        if self.cache and not os.path.exists(self.cache):
            logger.info(f"Created LCA cache dir: {self.cache}")
            os.makedirs(self.cache)
            self.tour_fp = os.path.join(self.cache, "tour.pkl.gz")
            self.level_fp = os.path.join(self.cache, "level.pkl.gz")
            self.occurrence_fp = os.path.join(self.cache, "occurrence.pkl.gz")
            self.sparse_fp = os.path.join(self.cache, "precomputed_lcas.pkl.gz")
        self.tour = None
        self.level = None
        self.occurrence = None
        self.sparse = None
        self.lca_prepared = False

    def prepare_tree(self):
        """Performs Eulerian tour of nodes.dmp taxids and constructs three data
        structures:

        1. tour : list of branches and leaves.
        2. level: list of distances from the root.
        3. occurrence: dict of occurences of the taxid respective to the root.

        Notes
        -----
            For more information on why we construct these three data structures see references below:

            * `Geeksforgeeks: Find LCA in Binary Tree using RMQ  https://www.geeksforgeeks.org/find-lca-in-binary-tree-using-rmq/`_
            * `Topcoder: Another easy solution in <O(N logN, O(logN)> https://www.topcoder.com/community/competitive-programming/tutorials/range-minimum-query-and-lowest-common-ancestor/#Another%20easy%20solution%20in%20O(N%20logN,%20O(logN)`_

        Returns
        -------
        NoneType
            sets internals to be used for LCA lookup

        """
        repickle = False
        # Creating a repickle toggle b/c if *any* of the file paths do not exist, we need to redo
        for fp in [self.tour_fp, self.level_fp, self.occurrence_fp]:
            if not os.path.exists(fp):
                repickle = True
        if not repickle:
            try:
                self.tour = unpickle(fpath=self.tour_fp)
                self.level = unpickle(fpath=self.level_fp)
                self.occurrence = unpickle(fpath=self.occurrence_fp)
                return
            except UnpicklingError as err:
                logger.error(
                    "Error unpickling `tour`,`level`, or `occurrence`. One of these files may be corrupted! Rebuilding..."
                )
                # This will now continue resulting in overwriting the corrupted files.
        if self.verbose:
            logger.debug("Preparing tree, level, occurrence for LCA/RMQ")
        taxids = {}
        parents = {}
        children = {}
        for taxid, info in self.nodes.items():
            if taxid == 1:
                # Skip root
                continue
            parent = info["parent"]
            parents.update({taxid: parent})
            taxids.update({taxid: 1})
            if parent in children:
                children[parent].add(taxid)
            else:
                children.update({parent: set([taxid])})
        # Start tour with root as (0,1)
        tour = [(0, 1)]
        direction = 1
        dist = 0
        level = [dist]
        while taxids:
            if direction > 0:
                parent = tour[-1][1]
                if parent not in children:
                    direction *= -1
                    continue
                child = children[parent].pop()
                new_node = (parent, child)
                tour.append(new_node)
                dist += 1
                level.append(dist)
                taxids.pop(child, None)
                # If the set is now empty, we need to delete the parent key in children
                if not children[parent]:
                    children.pop(parent, None)
                # Do nothing else
            elif direction < 0:
                child = tour[-1][1]
                if child not in parents:
                    direction *= -1
                    continue
                parent = parents[child]
                new_node = (child, parent)
                tour.append(new_node)
                dist -= 1
                level.append(dist)
                # Delete the child from the parents dictionary
                parents.pop(child, None)
                # If parent still has children, reverse direction
                if parent in children:
                    direction *= -1
        occurrence = {}
        for i, node in enumerate(tour):
            child = node[1]
            if child not in occurrence:
                occurrence[child] = i
        self.tour = tour
        self.level = level
        self.occurrence = occurrence
        if self.cache:
            make_pickle(obj=self.tour, outfpath=self.tour_fp)
            make_pickle(obj=self.level, outfpath=self.level_fp)
            make_pickle(obj=self.occurrence, outfpath=self.occurrence_fp)
        return

    def preprocess_minimums(self):
        """Preprocesses all possible LCAs.

        This constructs a sparse table to be used for LCA/Range Minimum Query
        using the `self.level` array associated with its respective eulerian `self.tour`.
        For more information on these data structures see :func:`~lca.LCA.prepare_tree`.

        Sparse table size:

        * n = number of elements in level list
        * rows range = (0 to n)
        * columns range = (0 to logn)

        Returns
        -------
        NoneType
            sets `self.sparse` internal to be used for LCA lookup.

        """
        if os.path.exists(self.sparse_fp) and os.path.getsize(self.sparse_fp):
            try:
                self.sparse = unpickle(fpath=self.sparse_fp)
                return
            except UnpicklingError as err:
                logger.error(
                    f"Error unpickling {self.sparse_fp}, this may be corrupted! Overwriting..."
                )
        logger.debug("Precomputing LCAs")
        # Instantiate an empty sparse array with dimensions from `self.level`
        nrows = len(self.level)
        ncols = int(np.floor(np.log2(nrows)) + 1)
        sparse_array = np.empty((nrows, ncols))
        sparse_array[:, 0] = [i for i in range(nrows)]
        # We start at 1th column because we have the indices in the 0th column from above
        for col in tqdm(
            range(1, ncols), disable=self.disable, desc="Precomputing LCAs", leave=False
        ):
            for row in range(0, nrows):
                # First we check that the exponent of the column does not exceed the rows
                if 2 ** col > nrows:
                    continue
                # Next check whether element at pos is within rows
                if row + (2 ** col) - 1 >= nrows:
                    sparse_array[row, col] = False
                    continue
                # We now have our range in terms of indices
                # Retrieve indices corresponding to unique LCA comparisons
                lower_index = sparse_array[row, (col - 1)]
                upper_index = sparse_array[(row + 2 ** (col - 1)), (col - 1)]
                # We need to cast ints here to conver numpy.float64 to ints for index accession.
                lower_index, upper_index = map(int, [lower_index, upper_index])
                # Access levels via comparison indices
                lower_min = self.level[lower_index]
                upper_min = self.level[upper_index]
                # Set the winning minimum between the range of indices
                if lower_min < upper_min:
                    sparse_array[row, col] = lower_index
                else:
                    sparse_array[row, col] = upper_index
        if self.cache:
            make_pickle(obj=sparse_array, outfpath=self.sparse_fp)
        self.sparse = sparse_array
        return

    def prepare_lca(self):
        """Prepare LCA internal data structures for :func:`~lca.LCA.lca`.

        e.g. self.tour, self.level, self.occurrence, self.sparse are all ready.

        Returns
        -------
        NoneType
            Prepares all LCA internals and if successful sets `self.lca_prepared` to True.

        """
        self.prepare_tree()
        self.preprocess_minimums()
        self.lca_prepared = True
        logger.debug(f"LCA data structures prepared")
        # tour, level, occurrence, sparse all ready
        return

    def lca(self, node1, node2):
        """Performs Range Minimum Query between 2 taxids.

        Parameters
        ----------
        node1 : int
            taxid

        node2 : int
            taxid

        Returns
        -------
        int
            LCA taxid

        Raises
        -------
        ValueError
            Provided taxid is not in the nodes.dmp tree.

        """
        if not self.lca_prepared:
            self.prepare_lca()
        if node1 == None and node2 == None:
            return 1
        if node1 not in self.occurrence:
            raise ValueError(f"{node1} not in tree")
        if node2 not in self.occurrence:
            raise ValueError(f"{node2} not in tree")
        if node1 == None:
            return node2
        if node2 == None:
            return node1
        if node1 == node2:
            return node1
        if self.occurrence[node1] < self.occurrence[node2]:
            low = self.occurrence[node1]
            high = self.occurrence[node2]
        else:
            low = self.occurrence[node2]
            high = self.occurrence[node1]
        # equipartition range b/w both nodes.
        cutoff_range = int(np.floor(np.log2(high - low + 1)))
        lower_index = self.sparse[low, cutoff_range]
        upper_index = self.sparse[(high - (2 ** cutoff_range) + 1), cutoff_range]
        lower_index, upper_index = map(int, [lower_index, upper_index])
        lower_range = self.level[lower_index]
        upper_range = self.level[upper_index]
        if lower_range <= upper_range:
            lca_range = lower_range
        else:
            lca_range = upper_range
        lca_node = self.tour[self.level.index(lca_range, low, high)]
        # (parent, child)
        return lca_node[1]

    def convert_sseqids_to_taxids(
        self,
        sseqids: Dict[str, Set[str]],
    ) -> Tuple[Dict[str, Set[int]], pd.DataFrame]:
        """
        Translates subject sequence ids to taxids from prot.accession2taxid.gz and dead_prot.accession2taxid.gz.

        .. note::

            If an accession number is no longer available in prot.accesssion2taxid.gz
            (either due to being suppressed, deprecated or removed by NCBI),
            then root taxid (1) is returned as the taxid for the corresponding sseqid.

        Parameters
        ----------
        sseqids : dict
            {qseqid: {sseqid, ...}, ...}

        Returns
        -------
        Tuple[Dict[str, Set[int]], pd.DataFrame]
            {qseqid: {taxid, taxid, ...}, ...}, index=range, cols=[qseqid, sseqid, raw_taxid, merged_taxid, cleaned_taxid]

        Raises
        -------
        FileNotFoundError
            prot.accession2taxid.gz database is required for sseqid to taxid conversion.

        """
        # We first get all unique sseqids that were retrieved from the blast output
        recovered_sseqids = set(
            chain.from_iterable([qseqid_sseqids for qseqid_sseqids in sseqids.values()])
        )

        # Check for sseqid in dead_prot.accession2taxid.gz in case an old db was used.
        # Any accessions not found in live prot.accession2taxid db *may* be here.
        # This *attempts* to prevent sseqids from being assigned root (root taxid=1)
        try:
            sseqids_to_taxids = self.search_prot_accessions(
                accessions=recovered_sseqids,
                db="dead",
                sseqids_to_taxids=None,
            )
            dead_sseqids_found = set(sseqids_to_taxids.keys())
        except FileNotFoundError as db_fpath:
            logger.warn(f"Skipping taxid conversion from {db_fpath}")
            sseqids_to_taxids = None
            dead_sseqids_found = set()

        # Now build the mapping from sseqid to taxid from the full/live accession2taxid dbs
        # Possibly overwriting any merged accessions to live accessions
        sseqids_to_taxids = self.search_prot_accessions(
            accessions=recovered_sseqids,
            db="full",
            sseqids_to_taxids=sseqids_to_taxids,
        )

        # Remove sseqids: Ignore any sseqids already found
        live_sseqids_found = set(sseqids_to_taxids.keys())
        live_sseqids_found -= dead_sseqids_found
        recovered_sseqids -= live_sseqids_found
        if recovered_sseqids:
            logger.warn(
                f"sseqids without corresponding taxid: {len(recovered_sseqids):,}"
            )
        root_taxid = 1
        taxids = {}
        sseqid_not_found = pd.NA
        sseqid_to_taxid_df = pd.DataFrame(
            [
                {
                    "qseqid": qseqid,
                    "sseqid": sseqid,
                    "raw_taxid": sseqids_to_taxids.get(sseqid, sseqid_not_found),
                }
                for qseqid, qseqid_sseqids in sseqids.items()
                for sseqid in qseqid_sseqids
            ]
        )
        # Update taxids if they are old.
        sseqid_to_taxid_df["merged_taxid"] = sseqid_to_taxid_df.raw_taxid.map(
            lambda tid: self.merged.get(tid, tid)
        )
        # If we still have missing taxids, we will set the sseqid value to the root taxid
        # fill missing taxids with root_taxid
        sseqid_to_taxid_df["cleaned_taxid"] = sseqid_to_taxid_df.merged_taxid.fillna(
            root_taxid
        )
        for qseqid, qseqid_sseqids in sseqids.items():
            # NOTE: we only want to retrieve the set of unique taxids (not a list) for LCA query
            qseqid_taxids = {
                sseqids_to_taxids.get(sseqid, root_taxid) for sseqid in qseqid_sseqids
            }
            taxids[qseqid] = qseqid_taxids
        return taxids, sseqid_to_taxid_df

    def reduce_taxids_to_lcas(
        self, taxids: Dict[str, Set[int]]
    ) -> Tuple[Dict[str, int], pd.DataFrame]:
        """Retrieves the lowest common ancestor for each set of taxids in of the taxids

        Parameters
        ----------
        taxids : dict
            {qseqid: {taxid, ...}, qseqid: {taxid, ...}, ...}

        Returns
        -------
        Tuple[Dict[str, int], pd.DataFrame]
            {qseqid: lca, qseqid: lca, ...}, pd.DataFrame(index=range, cols=[qseqid, taxids])

        """
        lca_taxids = {}
        missing_taxids = []
        root_taxid = 1
        logger.debug(f"Assigning LCAs to {len(taxids):,} ORFs")
        for qseqid, qseqid_taxids in taxids.items():
            # This will convert any deprecated taxids to their current taxids before reduction to LCA.
            # NOTE: we only want to retrieve the set of unique taxids (not a list) for LCA query
            qseqid_taxids = {self.merged.get(taxid, taxid) for taxid in qseqid_taxids}
            lca_taxid = False
            num_taxids = len(qseqid_taxids)
            while not lca_taxid:
                # Reduce all qseqid taxids to LCA
                if num_taxids >= 2:
                    try:
                        lca_taxid = functools.reduce(
                            lambda taxid1, taxid2: self.lca(node1=taxid1, node2=taxid2),
                            qseqid_taxids,
                        )
                    except ValueError as err:
                        logger.warn(
                            f"Missing taxids during LCA retrieval. Setting {qseqid} ({len(qseqid_taxids):,} ORF taxids) to root taxid."
                        )
                        lca_taxid = root_taxid
                        missing_taxids.append(
                            {
                                "qseqid": qseqid,
                                "taxids": ",".join(
                                    [str(taxid) for taxid in qseqid_taxids]
                                ),
                            }
                        )
                # If only one taxid was recovered... This is our LCA
                elif num_taxids == 1:
                    lca_taxid = qseqid_taxids.pop()
                else:
                    # Exception handling where input for qseqid contains no taxids e.g. num_taxids == 0:
                    lca_taxid = root_taxid
            lca_taxids.update({qseqid: lca_taxid})
        missing_df = pd.DataFrame(missing_taxids)
        return lca_taxids, missing_df

    def read_sseqid_to_taxid_table(
        self, sseqid_to_taxid_filepath: str
    ) -> Dict[str, Set[int]]:
        """Retrieve each qseqid's set of taxids from `sseqid_to_taxid_filepath` for reduction by LCA

        Parameters
        ----------
        sseqid_to_taxid_filepath : str
            Path to sseqid to taxid table with columns: qseqid, sseqid, raw_taxid, merged_taxid, clean_taxid

        Returns
        -------
        Dict[str, Set[int]]
            Dictionary keyed by qseqid containing sets of respective `clean` taxid
        """
        taxids = {}
        with open(sseqid_to_taxid_filepath) as fh:
            __ = fh.readline()
            for line in fh:
                qseqid, *__, clean_taxid = line.strip().split("\t")
                clean_taxid = int(clean_taxid)
                if qseqid in taxids:
                    taxids[qseqid].add(clean_taxid)
                else:
                    taxids[qseqid] = set([clean_taxid])
        return taxids

    def blast2lca(
        self,
        blast: str,
        out: str,
        sseqid_to_taxid_output: str = "",
        lca_reduction_log: str = "",
        force: bool = False,
    ) -> str:
        """Determine lowest common ancestor of provided amino-acid ORFs.

        Parameters
        ----------
        blast : str
            </path/to/diamond/outfmt6/blastp.tsv>.
        out : str
            </path/to/output/lca.tsv>.
        sseqid_to_taxid_output : str
            Path to write qseqids' sseqids and their taxid designations from NCBI databases
        force : bool, optional
            Force overwrite of existing `out`.

        Returns
        -------
        str
            `out` </path/to/output/lca.tsv>.

        """
        if os.path.exists(out) and os.path.getsize(out) and not force:
            logger.warning(f"FileAlreadyExists {out}")
            return out
        blast = os.path.realpath(blast)
        # If sseqid_to_taxid_output exists then we'll retrieve sseqid taxids from this...
        if (
            sseqid_to_taxid_output
            and os.path.exists(sseqid_to_taxid_output)
            and os.path.getsize(sseqid_to_taxid_output)
        ):
            logger.debug(f"Retrieving taxids from {sseqid_to_taxid_output}")
            taxids = self.read_sseqid_to_taxid_table(sseqid_to_taxid_output)
        elif not os.path.exists(blast) or not os.path.getsize(blast):
            raise FileNotFoundError(blast)
        elif os.path.exists(blast) and os.path.getsize(blast):
            logger.debug(f"Retrieving taxids from {blast}")
            sseqids = diamond.parse(results=blast, verbose=self.verbose)
            taxids, sseqid_to_taxid_df = self.convert_sseqids_to_taxids(sseqids)
            if sseqid_to_taxid_output:
                sseqid_to_taxid_df.to_csv(
                    sseqid_to_taxid_output, sep="\t", index=False, header=True
                )
                logger.info(
                    f"Wrote {sseqid_to_taxid_df.shape[0]:,} qseqids to {sseqid_to_taxid_output}"
                )
        # TODO: merge sseqid_to_taxid_output table with err_qseqids_df
        # Alluvial diagrams may be interesting here...
        lcas, err_qseqids_df = self.reduce_taxids_to_lcas(taxids)
        if lca_reduction_log:
            err_qseqids_df.to_csv(lca_reduction_log, sep="\t", index=False, header=True)
        written_lcas = self.write_lcas(lcas=lcas, out=out)
        return written_lcas

    def write_lcas(self, lcas: Dict[str, int], out: str) -> str:
        """Write `lcas` to tab-delimited file: `out`.

        Ordered columns are:

            * qseqid : query seqid
            * name : LCA name
            * rank : LCA rank
            * lca : LCA taxid

        Parameters
        ----------
        lcas : dict
            {qseqid:lca_taxid, qseqid:lca_taxid, ...}
        out : str
            </path/to/output/file.tsv>

        Returns
        -------
        str
            `out`

        """
        lines = "qseqid\tname\trank\tlca\n"
        fh = open(out, "w")
        nlines = 0
        buffer_size = 10000
        for qseqid, taxid in lcas.items():
            if nlines >= buffer_size:
                fh.write(lines)
                lines = ""
                nlines = 0
            lines += f"{qseqid}\t{self.name(taxid)}\t{self.rank(taxid)}\t{taxid}\n"
            nlines += 1
        fh.write(lines)
        fh.close()
        logger.info(f"Assigned LCA to {len(lcas):,} ORFs: {out}")
        return out

    def parse(
        self, lca_fpath: str, orfs_fpath: str = None
    ) -> Dict[str, Dict[str, Dict[int, int]]]:
        """Retrieve and construct contig dictionary from provided `lca_fpath`.

        Parameters
        ----------
        lca_fpath : str
            </path/to/lcas.tsv>
            tab-delimited ordered columns: qseqid, name, rank, lca_taxid

        orfs_fpath : str, optional (required if using prodigal version <2.6)
            </path/to/prodigal/called/orfs.fasta>
            Note: These ORFs should correspond to the ORFs provided in the BLAST table.

        Returns
        -------
        dict
            {contig:{rank:{taxid:counts, ...}, rank:{...}, ...}, ...}

        Raises
        -------
        FileNotFoundError
            `lca_fpath` does not exist.

        FileNotFoundError
            `orfs_fpath` does not exist.

        ValueError
            If prodigal version is under 2.6, `orfs_fpath` is a required input.

        """
        logger.debug(f"Parsing LCA table: {lca_fpath}")
        if not os.path.exists(lca_fpath):
            raise FileNotFoundError(lca_fpath)
        if orfs_fpath and not os.path.exists(orfs_fpath):
            raise FileNotFoundError(orfs_fpath)

        version = get_versions("prodigal")
        if version.count(".") >= 2:
            version = float(".".join(version.split(".")[:2]))
        else:
            version = float(version)
        if version < 2.6 and not orfs_fpath:
            raise ValueError("Prodigal version under 2.6 requires orfs_fpath input!")
        # Create a contig translation dictionary with or without ORFs
        if orfs_fpath:
            contigs_from_orfs = prodigal.contigs_from_headers(orfs_fpath)
        else:
            df = pd.read_csv(lca_fpath, sep="\t", usecols=["qseqid"])
            df["contig"] = df.qseqid.map(lambda orf: orf.rsplit("_", 1)[0])
            contigs_from_orfs = df.set_index("qseqid")["contig"].to_dict()

        fname = os.path.basename(lca_fpath)
        n_lines = file_length(lca_fpath) if self.verbose else None
        lca_hits = {}
        with open(lca_fpath) as fh:
            __ = fh.readline()  # header
            for line in tqdm(
                fh,
                total=n_lines,
                disable=self.disable,
                desc=f"Parsing {fname}",
                leave=False,
            ):
                # orf, name, rank, taxid
                orf_id, __, rank, taxid = line.strip().split("\t")
                taxid = int(taxid)
                contig = contigs_from_orfs.get(orf_id)
                if taxid != 1:
                    while rank not in set(NCBI.CANONICAL_RANKS):
                        taxid = self.parent(taxid)
                        rank = self.rank(taxid)
                if contig not in lca_hits:
                    lca_hits.update({contig: {rank: {taxid: 1}}})
                elif rank not in lca_hits[contig]:
                    lca_hits[contig].update({rank: {taxid: 1}})
                elif taxid not in lca_hits[contig][rank]:
                    lca_hits[contig][rank].update({taxid: 1})
                else:
                    lca_hits[contig][rank][taxid] += 1
        return lca_hits


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Script to determine Lowest Common Ancestor",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--blast",
        help="Path to BLAST results table respective to `orfs`. "
        "(Note: The table provided must be in outfmt=6)",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--dbdir",
        help="Path to NCBI databases directory.",
        metavar="dirpath",
        default=NCBI_DIR,
    )
    parser.add_argument(
        "--lca-output",
        help="Path to write LCA results.",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--sseqid2taxid-output",
        help="Path to write qseqids sseqids to taxids translations table",
        required=False,
        metavar="filepath",
    )
    parser.add_argument(
        "--lca-error-taxids",
        help="Path to write table of blast table qseqids that were assigned root due to a missing taxid",
        required=False,
        metavar="filepath",
    )
    parser.add_argument(
        "--verbose",
        help="Add verbosity to logging stream.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--force",
        help="Force overwrite if results already exist.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--cache", help="Path to cache pickled LCA database objects.", metavar="dirpath"
    )
    parser.add_argument(
        "--only-prepare-cache",
        help="Only prepare the LCA database objects and write to provided --cache parameter",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--force-cache-overwrite",
        help="Force overwrite if results already exist.",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    lca = LCA(dbdir=args.dbdir, verbose=args.verbose, cache=args.cache)

    # The pkl objects will be recomputed if the respective process cache file does not exist
    if args.force_cache_overwrite:
        logger.info("Overwriting cache")
        if os.path.exists(lca.tour_fp):
            os.remove(lca.tour_fp)
        if os.path.exists(lca.sparse_fp):
            os.remove(lca.sparse_fp)

    if args.only_prepare_cache:
        lca.prepare_lca()
        logger.info("Cache prepared")
        return

    lca.blast2lca(
        blast=args.blast,
        out=args.lca_output,
        sseqid_to_taxid_output=args.sseqid2taxid_output,
        lca_reduction_log=args.lca_error_taxids,
        force=args.force,
    )


if __name__ == "__main__":
    main()
