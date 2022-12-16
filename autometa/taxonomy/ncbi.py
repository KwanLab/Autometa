#!/usr/bin/env python
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

File containing definition of the NCBI class and containing functions useful for handling NCBI taxonomy databases
"""


import gzip
import logging
import os
import string
import sys

from typing import Dict, Set, Tuple
from itertools import chain

import pandas as pd

from tqdm import tqdm

from autometa.common.utilities import file_length
from autometa.common.exceptions import DatabaseOutOfSyncError
from autometa.config.utilities import DEFAULT_CONFIG
from autometa.taxonomy.database import TaxonomyDatabase


logger = logging.getLogger(__name__)

NCBI_DIR = DEFAULT_CONFIG.get("databases", "ncbi")
# For cases where autometa has not been configured, attempt to find ncbi directory via source
NCBI_DIR = NCBI_DIR if not "None" in NCBI_DIR else NCBI_DIR.replace("None", ".")


class NCBI(TaxonomyDatabase):
    """Taxonomy utilities for NCBI databases."""

    def __init__(self, dbdir, verbose=False):
        """
        Instantiates the NCBI class

        Parameters
        ----------

        dbdir : str
            Path to the database directory

        verbose : bool, optional
            log progress to terminal, by default False

        """
        self.dbdir = dbdir
        self.verbose = verbose
        self.disable = not self.verbose
        self.names_fpath = os.path.join(self.dbdir, "names.dmp")
        self.nodes_fpath = os.path.join(self.dbdir, "nodes.dmp")
        self.merged_fpath = os.path.join(self.dbdir, "merged.dmp")
        self.delnodes_fpath = os.path.join(self.dbdir, "delnodes.dmp")
        # Set prot.accession2taxid filepath
        self.accession2taxid_fpath = os.path.join(self.dbdir, "prot.accession2taxid")
        acc2taxid_gz = ".".join([self.accession2taxid_fpath, "gz"])
        if not os.path.exists(self.accession2taxid_fpath) and os.path.exists(
            acc2taxid_gz
        ):
            self.accession2taxid_fpath = acc2taxid_gz
        # Set prot.accession2taxid.FULL.gz filepath
        self.accession2taxidfull_fpath = os.path.join(
            self.dbdir, "prot.accession2taxid.FULL"
        )
        acc2taxid_gz = ".".join([self.accession2taxidfull_fpath, "gz"])
        if not os.path.exists(self.accession2taxidfull_fpath) and os.path.exists(
            acc2taxid_gz
        ):
            self.accession2taxidfull_fpath = acc2taxid_gz
        # Set dead_prot.accession2taxid filepath
        self.dead_accession2taxid_fpath = os.path.join(
            self.dbdir, "dead_prot.accession2taxid"
        )
        acc2taxid_gz = ".".join([self.dead_accession2taxid_fpath, "gz"])
        if not os.path.exists(self.dead_accession2taxid_fpath) and os.path.exists(
            acc2taxid_gz
        ):
            self.dead_accession2taxid_fpath = acc2taxid_gz
        # Check formatting for nr database
        self.nr_fpath = os.path.join(self.dbdir, "nr.gz")
        nr_bname = os.path.splitext(os.path.basename(self.nr_fpath))[0]
        nr_dmnd_fname = ".".join([nr_bname, "dmnd"])
        nr_dmnd_fpath = os.path.join(self.dbdir, nr_dmnd_fname)
        if os.path.exists(nr_dmnd_fpath):
            self.nr_fpath = nr_dmnd_fpath
        if "dmnd" not in os.path.basename(self.nr_fpath):
            # This check should probably be in the dependencies/databases file...
            logger.warning(
                f"DatabaseWarning: {self.nr_fpath} needs to be formatted for diamond!"
            )
        # Setup data structures from nodes.dmp, names.dmp and merged.dmp
        self.nodes = self.parse_nodes()
        self.names = self.parse_names()
        self.merged = self.parse_merged()
        self.delnodes = self.parse_delnodes()

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
        return self.dbdir

    def parse_names(self) -> Dict[int, str]:
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

    def parse_nodes(self) -> Dict[int, str]:
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
            parent, child = [int(node) for node in [parent, child]]
            rank = rank.lower()
            nodes.update({child: {"parent": parent, "rank": rank}})
        fh.close()
        if self.verbose:
            logger.debug("nodes loaded")
        return nodes

    def parse_merged(self) -> Dict[int, int]:
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

    def parse_delnodes(self) -> Set[int]:
        """
        Parse the delnodes.dmp database
        Note: This is performed when a new NCBI class instance is constructed

        Returns
        -------
        set
            {taxid, ...}
        """
        if self.verbose:
            logger.debug(f"Processing delnodes from {self.delnodes_fpath}")
        fh = open(self.delnodes_fpath)
        delnodes = set()
        for line in tqdm(
            fh, disable=self.disable, desc="parsing delnodes", leave=False
        ):
            del_taxid = int(line.strip("\t|\n"))
            delnodes.add(del_taxid)
        fh.close()
        if self.verbose:
            logger.debug("delnodes loaded")
        return delnodes

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

    def search_prot_accessions(
        self,
        accessions: set,
        sseqids_to_taxids: Dict[str, int] = None,
        db: str = "live",
    ) -> Dict[str, int]:
        """Search prot.accession2taxid.gz and dead_prot.accession2taxid.gz

        Parameters
        ----------
        accessions : set
            Set of subject sequence ids retrieved from diamond blastp search (sseqids)

        sseqids_to_taxids : Dict[str, int], optional
            Dictionary containing sseqids converted to taxids

        db : str, optional
            selection of one of the prot accession to taxid databases from NCBI. Choices are live, dead, full

            * live: prot.accession2taxid.gz
            * full: prot.accession2taxid.FULL.gz
            * dead: dead_prot.accession2taxid.gz

        Returns
        -------
        Dict[str, int]
            Dictionary containing sseqids converted to taxids
        """
        if not sseqids_to_taxids:
            sseqids_to_taxids = {}
        if not isinstance(db, str):
            raise ValueError(f"db must be a string! Type Given: {type(db)}")
        db = db.lower()
        choices = {
            "live": self.accession2taxid_fpath,
            "dead": self.dead_accession2taxid_fpath,
            "full": self.accession2taxidfull_fpath,
        }
        if db not in choices:
            raise ValueError(f"db must be one of live, full or dead. Given: {db}")
        fpath = choices.get(db)

        # Revert to accession2taxid if FULL is not present
        if db == "full" and (not os.path.exists(fpath) or not os.path.getsize(fpath)):
            logger.warn(
                "prot.accession2taxid.FULL.gz was not found. Reverting to prot.accession2taxid.gz"
            )
            logger.warn(
                "To achieve greater resolution of your metagenome taxonomy, considering downloading the prot.accession2taxid.FULL.gz database file"
            )
            fpath = choices.get("live")
            db = "live"

        if not os.path.exists(fpath) or not os.path.getsize(fpath):
            raise FileNotFoundError(fpath)

        # "rt" open the database in text mode instead of binary to be handled like a text file
        fh = gzip.open(fpath, "rt") if fpath.endswith(".gz") else open(fpath)
        filename = os.path.basename(fpath)
        # skip the header line
        __ = fh.readline()
        logger.debug(
            f"Searching for {len(accessions):,} accessions in {filename}. This may take a while..."
        )
        n_lines = file_length(fpath, approximate=True) if self.verbose else None
        desc = f"Parsing {filename}"
        converted_sseqid_count = 0
        for line in tqdm(
            fh, disable=self.disable, desc=desc, total=n_lines, leave=False
        ):
            if db == "full":
                # FULL format is accession.version\ttaxid\n
                acc_num = None  # Just in case
                acc_ver, taxid = line.strip().split("\t")
            else:
                # dead and live formats are accession\taccession.version\ttaxid\tgi\n
                acc_num, acc_ver, taxid, _ = line.strip().split("\t")

            taxid = int(taxid)
            if acc_ver in accessions:
                sseqids_to_taxids[acc_ver] = taxid
                converted_sseqid_count += 1

            # So prog will not have to search through the accessions set
            if db == "full":
                continue

            # Search for base accession if using live or dead accession2taxid databases
            if acc_num in accessions:
                sseqids_to_taxids[acc_num] = taxid
                converted_sseqid_count += 1

        fh.close()
        logger.debug(f"sseqids converted from {filename}: {converted_sseqid_count:,}")
        return sseqids_to_taxids

    def convert_accessions_to_taxids(
        self, accessions: set
    ) -> Tuple[Dict[str, Set[int]], pd.DataFrame]:
        # We first get all unique accessions that were retrieved from the blast output
        recovered_accessions = set(
            chain.from_iterable(
                [qseqid_sseqids for qseqid_sseqids in accessions.values()]
            )
        )
        # Check for sseqid in dead_prot.accession2taxid.gz in case an old db was used.
        # Any accessions not found in live prot.accession2taxid db *may* be here.
        # This *attempts* to prevent accessions from being assigned root (root taxid=1)
        try:
            sseqids_to_taxids = self.search_prot_accessions(
                accessions=recovered_accessions,
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
            accessions=recovered_accessions,
            db="full",
            sseqids_to_taxids=sseqids_to_taxids,
        )

        # Remove accessions: Ignore any accessions already found
        live_sseqids_found = set(sseqids_to_taxids.keys())
        live_sseqids_found -= dead_sseqids_found
        recovered_accessions -= live_sseqids_found
        if recovered_accessions:
            logger.warn(
                f"accessions without corresponding taxid: {len(recovered_accessions):,}"
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
                for qseqid, qseqid_sseqids in accessions.items()
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
        root_taxid = 1
        for qseqid, qseqid_sseqids in accessions.items():
            # NOTE: we only want to retrieve the set of unique taxids (not a list) for LCA query
            qseqid_taxids = {
                sseqids_to_taxids.get(sseqid, root_taxid) for sseqid in qseqid_sseqids
            }
            taxids[qseqid] = qseqid_taxids
        return taxids, sseqid_to_taxid_df


if __name__ == "__main__":
    print("file containing Autometa NCBI utilities class")
    import sys

    sys.exit(1)
