#!/usr/bin/env python
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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

This script contains the LCA class containing methods to determine the
Lowest Common Ancestor given a tab-delimited BLAST table, fasta file, or
iterable of SeqRecords.

Note: LCA will assume the BLAST results table is in output format 6.
"""


import functools
import gzip
import logging
import os

from typing import Dict
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

    def __init__(self, dbdir: str, verbose: bool = False, outdir: str = None):
        super().__init__(dbdir, verbose=verbose)
        self.verbose = verbose
        self.dbdir = dbdir
        self.outdir = outdir
        self.disable = False if verbose else True
        self.tour_fp = os.path.join(self.outdir, "tour.pkl.gz")
        self.tour = None
        self.level_fp = os.path.join(self.outdir, "level.pkl.gz")
        self.level = None
        self.occurrence_fp = os.path.join(self.outdir, "occurrence.pkl.gz")
        self.occurrence = None
        self.sparse_fp = os.path.join(self.outdir, "precomputed_lcas.pkl.gz")
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
        if self.outdir:
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
        if self.outdir:
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

    def search_prot_accessions(
        self,
        accessions: set,
        sseqids_to_taxids: Dict[str, int] = None,
        live: bool = True,
    ) -> Dict[str, int]:
        """Search prot.accession2taxid.gz and dead_prot.accession2taxid.gz

        Parameters
        ----------
        accessions : set
            Set of subject sequence ids retrieved from diamond blastp search (sseqids)

        sseqids_to_taxids : Dict[str, int], optional
            Dictionary containing sseqids converted to taxids

        live : bool, optional
            Whether to use prot.accession2taxid.gz (live=True) or dead_prot.accession2taxid.gz (live=False)

        Returns
        -------
        Dict[str, int]
            Dictionary containing sseqids converted to taxids
        """
        if not sseqids_to_taxids:
            sseqids_to_taxids = {}
        if live:
            # "rt" open the database in text mode instead of binary to be handled like a text file
            fh = (
                gzip.open(self.accession2taxid_fpath, "rt")
                if self.accession2taxid_fpath.endswith(".gz")
                else open(self.accession2taxid_fpath)
            )
            fpath = self.accession2taxid_fpath
        else:
            fh = (
                gzip.open(self.dead_accession2taxid_fpath, "rt")
                if self.dead_accession2taxid_fpath.endswith(".gz")
                else open(self.dead_accession2taxid_fpath)
            )
            fpath = self.dead_accession2taxid_fpath
        filename = os.path.basename(fpath)
        # skip the header line
        __ = fh.readline()
        logger.debug(
            f"Searching for {len(accessions):,} accessions in {filename}. This may take a while..."
        )
        n_lines = file_length(fpath, approximate=True) if self.verbose else None
        desc = f"Parsing {filename}"
        sseqids_to_taxids = {}
        converted_sseqid_count = 0
        for line in tqdm(
            fh, disable=self.disable, desc=desc, total=n_lines, leave=False
        ):
            acc_num, acc_ver, taxid, _ = line.split("\t")
            taxid = int(taxid)
            if acc_num in accessions:
                sseqids_to_taxids[acc_num] = taxid
                converted_sseqid_count += 1
            if acc_ver in accessions:
                sseqids_to_taxids[acc_ver] = taxid
                converted_sseqid_count += 1
        fh.close()
        logger.debug(f"sseqids converted from {filename}: {converted_sseqid_count:,}")
        return sseqids_to_taxids

    def convert_sseqids_to_taxids(self, sseqids: Dict[str, str]) -> Dict[str, int]:
        """
        Translates subject sequence ids to taxids from prot.accession2taxid.gz and dead_prot.accession2taxid.gz.

        Note
        ----
        If an accession number is no longer available in prot.accesssion2taxid.gz
        (either due to being suppressed, deprecated or removed by NCBI),
        then root taxid (1) is returned as the taxid for the corresponding sseqid.

        Parameters
        ----------
        sseqids : dict
            {qseqid: {sseqid, ...}, ...}

        Returns
        -------
        dict
            {qseqid: {taxid, taxid, ...}, ...}

        Raises
        -------
        FileNotFoundError
            prot.accession2taxid.gz database is required for sseqid to taxid conversion.

        """
        # We first get all unique sseqids that were retrieved from the blast output
        recovered_sseqids = set(
            chain.from_iterable([qseqid_sseqids for qseqid_sseqids in sseqids.values()])
        )

        # Now build the mapping from sseqid to taxid from the live accession2taxid db
        sseqids_to_taxids = self.search_prot_accessions(
            accessions=recovered_sseqids,
            live=True,
            sseqids_to_taxids=None,
        )

        # Remove sseqids: Ignore any sseqids already found
        live_sseqids_found = set(sseqids_to_taxids.keys())
        recovered_sseqids -= live_sseqids_found
        # Next check for sseqid in dead_prot.accession2taxid.gz in case an old db was used.
        # Any accessions not found in live prot.accession2taxid db *may* be here.
        # This *attempts* to prevent sseqids from being assigned root (root taxid=1)
        sseqids_to_taxids = self.search_prot_accessions(
            accessions=recovered_sseqids,
            live=False,
            sseqids_to_taxids=sseqids_to_taxids,
        )
        dead_sseqids_found = set(sseqids_to_taxids.keys())
        dead_sseqids_found -= live_sseqids_found
        recovered_sseqids -= dead_sseqids_found
        root_taxid = 1
        taxids = {}
        for qseqid, qseqid_sseqids in recovered_sseqids.items():
            # NOTE: we only want to retrieve the set of unique taxids (not a list) for LCA query
            qseqid_taxids = {
                sseqids_to_taxids.get(sseqid, root_taxid) for sseqid in qseqid_sseqids
            }
            taxids[qseqid] = qseqid_taxids
        return taxids

    def reduce_taxids_to_lcas(self, taxids: Dict[str, set(int)]) -> Dict[str, int]:
        """Retrieves the lowest common ancestor for each set of taxids in of the taxids

        Parameters
        ----------
        taxids : dict
            {qseqid: {taxid, ...}, qseqid: {taxid, ...}, ...}

        Returns
        -------
        dict
            {qseqid: lca, qseqid: lca, ...}

        """
        lca_taxids = {}
        root_taxid = 1
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
                    except ValueError:
                        logger.error(
                            f"Missing taxids during LCA retrieval: {qseqid_taxids}. Setting {qseqid} to root taxid"
                        )
                        lca_taxid = root_taxid
                # If only one taxid was recovered... This is our LCA
                elif num_taxids == 1:
                    lca_taxid = qseqid_taxids.pop()
                else:
                    # Exception handling where input for qseqid contains no taxids e.g. num_taxids == 0:
                    lca_taxid = root_taxid
            lca_taxids.update({qseqid: lca_taxid})
        return lca_taxids

    def blast2lca(self, blast: str, out: str, force: bool = False) -> str:
        """Determine lowest common ancestor of provided amino-acid ORFs.

        Parameters
        ----------
        blast : str
            </path/to/diamond/outfmt6/blastp.tsv>.
        out : str
            </path/to/output/lca.tsv>.
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
        if not os.path.exists(blast) or not os.path.getsize(blast):
            raise FileNotFoundError(blast)
        blast = os.path.realpath(blast)
        sseqids = diamond.parse(results=blast, verbose=self.verbose)
        taxids = self.convert_sseqids_to_taxids(sseqids)
        lcas = self.reduce_taxids_to_lcas(taxids)
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

    parser = argparse.ArgumentParser(
        description="Script to determine Lowest Common Ancestor",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--output", help="Path to write LCA results.", required=True)
    parser.add_argument(
        "--blast",
        help="Path to BLAST results table respective to `orfs`. "
        "(Note: The table provided must be in outfmt=6)",
        required=True,
    )
    parser.add_argument(
        "--dbdir", help="Path to NCBI databases directory.", default=NCBI_DIR
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
    args = parser.parse_args()

    outdir = os.path.dirname(os.path.realpath(args.output))
    lca = LCA(dbdir=args.dbdir, verbose=args.verbose, outdir=outdir)

    lca.blast2lca(
        blast=args.blast,
        out=args.output,
        force=args.force,
    )


if __name__ == "__main__":
    main()
