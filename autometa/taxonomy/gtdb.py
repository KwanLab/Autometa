"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

File containing definition of the GTDB class and containing functions useful for handling GTDB taxonomy databases
"""


import shutil
import gzip
import logging
import os
import re
import string
import tarfile
import glob
from pathlib import Path

from typing import Dict, List, Set, Tuple
from itertools import chain
from tqdm import tqdm
from typing import Dict

import pandas as pd
import multiprocessing as mp

from autometa.common.utilities import file_length, is_gz_file
from autometa.common.external import diamond
from autometa.taxonomy.database import TaxonomyDatabase
from autometa.common.exceptions import DatabaseOutOfSyncError


logger = logging.getLogger(__name__)


def create_gtdb_db(reps_faa: str, dbdir: str) -> str:
    """
    Generate a combined faa file to create the GTDB-t database.

    Parameters
    ----------
    reps_faa : str
        Directory having faa file of all representative genomes. Can be tarballed.
    dbdir : str
        Path to output directory.

    Returns
    -------
    str
        Path to combined faa file. This can be used to make a diamond database.
    """

    if reps_faa.endswith(".tar.gz"):
        logger.debug(
            f"Extracting tarball containing GTDB ref genome animo acid data sequences to: {dbdir}/protein_faa_reps"
        )
        tar = tarfile.open(reps_faa)
        tar.extractall(path=dbdir)
        tar.close()
        logger.debug("Extraction done.")
        reps_faa = dbdir
 
    genome_protein_faa_filepaths = glob.glob(
        os.path.join(reps_faa, "**", "*_protein.faa.gz"), recursive=True
    )

    faa_index: Dict[str, str] = {}
    for genome_protein_faa_filepath in genome_protein_faa_filepaths:
        # Regex to get the genome accession from the file path
        genome_acc_search = re.search(
            r"GCA_\d+.\d?|GCF_\d+.\d?", genome_protein_faa_filepath
        )
        if genome_acc_search:
            faa_index[genome_protein_faa_filepath] = genome_acc_search.group()

    # Create dbdir if it doesn't exist
    if not os.path.isdir(dbdir):
        os.makedirs(dbdir)

    logger.debug(f"Merging {len(faa_index):,} faa files.")
    combined_faa = os.path.join(dbdir, "gtdb.faa") 
    with open(combined_faa, "w") as f_out:
        for faa_file, acc in faa_index.items():
            with gzip.open(faa_file, "rb") as f_in:
                for line in f_in:
                    # print(line.decode('utf-8'))
                    if line.decode('utf-8').startswith(">"):
                        seqheader = line.decode('utf-8').lstrip(">")
                        line = f"\n>{acc} {seqheader}"
                        f_out.write(line)
                    else:    
                        f_out.write(line.decode('utf-8'))
    logger.debug(f"Combined GTDB faa file written to {combined_faa}") 
    return combined_faa 


class GTDB(TaxonomyDatabase):
    """Taxonomy utilities for GTDB databases."""

    def __init__(self, dbdir: str, verbose: bool = True):
        """
        Instantiates the GTDB class

        """
        self.dbdir = dbdir
        self.verbose = verbose
        self.disable = not self.verbose
        self.dmnd_db = os.path.join(self.dbdir, "gtdb.dmnd")
        self.accession2taxid = os.path.join(self.dbdir, "taxid.map")
        self.nodes_fpath = os.path.join(self.dbdir, "nodes.dmp")
        self.names_fpath = os.path.join(self.dbdir, "names.dmp")
        self.merged_fpath = os.path.join(self.dbdir, "merged.dmp")
        self.delnodes_fpath = os.path.join(self.dbdir, "delnodes.dmp")
        self.verify_databases()
        self.names = self.parse_names()
        self.nodes = self.parse_nodes()
        self.merged = self.parse_merged()
        self.delnodes = self.parse_delnodes()

    def verify_databases(self):
        """
        Verify if the required databases are present.

        Raises
        ------
        FileNotFoundError
            One or more of the required database were not found.
        """
        for filepath in [
            self.names_fpath,
            self.nodes_fpath,
            self.dmnd_db,
            self.accession2taxid,
        ]:
            if not os.path.exists(filepath):
                raise FileNotFoundError(filepath)

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

    def search_genome_accessions(self, accessions: set) -> Dict[str, int]:
        """
        Search taxid.map file

        Parameters
        ----------
        accessions : set
            Set of subject sequence ids retrieved from diamond blastp search (sseqids)

        Returns
        -------
        Dict[str, int]
            Dictionary containing sseqids converted to taxids
        """
        sseqids_to_taxids = {}
        # "rt" open the database in text mode instead of binary to be handled like a text file
        fh = (
            gzip.open(self.accession2taxid, "rt")
            if is_gz_file(self.accession2taxid)
            else open(self.accession2taxid)
        )
        # skip the header line
        __ = fh.readline()
        n_lines = file_length(self.accession2taxid)
        converted_sseqid_count = 0
        logger.debug(f"parsing {self.accession2taxid}...")
        for line in fh:
            acc_ver, taxid = line.strip().split("\t")
            taxid = int(taxid)
            if acc_ver in accessions:
                sseqids_to_taxids[acc_ver] = taxid
                converted_sseqid_count += 1
        fh.close()
        logger.debug(f"sseqids converted: {converted_sseqid_count:,}")
        return sseqids_to_taxids

    def parse_nodes(self) -> Dict[int, str]:
        """
        Parse the `nodes.dmp` database.
        Note: This is performed when a new GTDB class instance is constructed

        Returns
        -------
        dict
            {child_taxid:{'parent':parent_taxid,'rank':rank}, ...}
        """
        fh = open(self.nodes_fpath)
        __ = fh.readline()  # root line
        nodes = {1: {"parent": 1, "rank": "root"}}
        for line in tqdm(fh, desc="parsing nodes", leave=False):
            child, parent, rank = line.split("\t|\t")[:3]
            parent, child = [int(node) for node in [parent, child]]
            rank = rank.lower()
            nodes.update({child: {"parent": parent, "rank": rank}})
        fh.close()
        return nodes

    def parse_names(self) -> Dict[int, str]:
        """
        Parses through names.dmp database and loads taxids with scientific names

        Returns
        -------
        dict
            {taxid:name, ...}
        """
        names = {}
        fh = open(self.names_fpath)
        for line in tqdm(fh, desc="parsing names", leave=False):
            taxid, name, __, classification = line.strip("\t|\n").split("\t|\t")[:4]
            taxid = int(taxid)
            name = name.lower()
            # Only add scientific name entries
            if classification == "scientific name":
                names.update({taxid: name})
        fh.close()
        return names

    def parse_merged(self) -> Dict[int, int]:
        """
        Parse the merged.dmp database

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

    def convert_accessions_to_taxids(
        self, accessions: Dict[str, Set[str]]
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
        recovered_accessions = set(
            chain.from_iterable(
                [qseqid_sseqids for qseqid_sseqids in accessions.values()]
            )
        )
        sseqids_to_taxids = self.search_genome_accessions(recovered_accessions)
        sseqid_not_found = pd.NA
        accession_to_taxid_df = pd.DataFrame(
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

        taxids = {}
        root_taxid = 1
        for qseqid, qseqid_sseqids in accessions.items():
            qseqid_taxids = {
                sseqids_to_taxids.get(sseqid, root_taxid) for sseqid in qseqid_sseqids
            }
            taxids[qseqid] = qseqid_taxids
        return taxids, accession_to_taxid_df


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--reps-faa",
        help="Path to directory containing GTDB ref genome animo acid data sequences. Can be tarballed.",
        required=True,
    )
    parser.add_argument(
        "--dbdir", help="Path to output GTDB database directory", required=True
    )
    parser.add_argument(
        "--cpus",
        help="Number of cpus to use for diamond-formatting GTDB database",
        default=mp.cpu_count(),
    )

    args = parser.parse_args()

    gtdb_combined = create_gtdb_db(reps_faa=args.reps_faa, dbdir=args.dbdir)
    diamond.makedatabase(
        fasta=gtdb_combined,
        database=gtdb_combined.replace(".faa", ".dmnd"),
        cpus=args.cpus,
    )


if __name__ == "__main__":
    # python -m autometa.taxonomy.gtdb --dbdir </path/to/gtdb>
    main()
