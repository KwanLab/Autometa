"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

File containing definition of the GTDB class and containing functions useful for handling GTDB taxonomy databases
"""

import argparse
import gzip
import logging
import os
import subprocess
import shutil

from typing import Dict, Iterable, List, Set, Tuple
from itertools import chain

import pandas as pd

from tqdm import tqdm

from autometa.common.utilities import file_length
from autometa.common.exceptions import DatabaseOutOfSyncError
from autometa.config.utilities import DEFAULT_CONFIG
from autometa.taxonomy.database import TaxonomyDatabase


logger = logging.getLogger(__name__)


def is_gz_file(filepath) -> bool:
    # https://stackoverflow.com/a/47080739
    with open(filepath, "rb") as test_f:
        return test_f.read(2) == b"\x1f\x8b"


def make_taxdump_files(taxa_files: List[str] = None, dbdir: str = None) -> None:
    cmd = ["gtdb_to_taxdump.py"]
    for taxa_file in taxa_files:
        cmd.append(taxa_file)
    if os.path.exists(dbdir):
        shutil.rmtree(dbdir)
    os.makedirs(dbdir)
    cmd.extend(["--outdir", dbdir])
    logger.debug(f'Running gtdb_to_taxdump.py: {" ".join(cmd)}')
    taxa_stdOut = os.path.join(dbdir, "taxID_info.tsv")
    fh = open(taxa_stdOut, "w")
    subprocess.run(cmd, stdout=fh, check=True)
    fh.close()


def make_gtdb_db(
    faa_tarball: str,
    names_fpath: str,
    nodes_fpath: str,
    dbdir: str,
    tmpdir: str = None,
    keep_temp: bool = False,
):
    dbdir = os.path.join(dbdir, "gtdb_to_diamond")
    os.makedirs(dbdir, exist_ok=True)
    cmd = [
        "gtdb_to_diamond.py",
        faa_tarball,
        names_fpath,
        nodes_fpath,
        "--outdir",
        dbdir,
    ]
    if tmpdir:
        cmd.extend(["--tmpdir", tmpdir])
    if keep_temp:
        cmd.append("--keep-temp")
    logger.debug(f'Running gtdb_to_diamond.py: {" ".join(cmd)}')
    subprocess.run(cmd, check=True)


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
        self.accession2taxid = os.path.join(self.dbdir, "accession2taxid.tsv")
        self.nodes_fpath = os.path.join(self.dbdir, "nodes.dmp")
        self.names_fpath = os.path.join(self.dbdir, "names.dmp")
        self.verify_databases()
        self.prepare_databases()

    def verify_databases(self):
        for filepath in [
            self.names_fpath,
            self.nodes_fpath,
            self.dmnd_db,
            self.accession2taxid,
        ]:
            if not os.path.exists(filepath):
                raise FileNotFoundError(filepath)

    def prepare_databases(self) -> None:
        self.names = self.parse_names()
        self.nodes = self.parse_nodes()

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
            acc_num, acc_ver, taxid = line.strip().split("\t")
            taxid = int(taxid)
            if acc_ver in accessions:
                sseqids_to_taxids[acc_ver] = taxid
                converted_sseqid_count += 1
        fh.close()
        logger.debug(f"sseqids converted: {converted_sseqid_count:,}")
        return sseqids_to_taxids

    def parse_nodes(self) -> Dict[int, str]:
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

    def convert_accessions_to_taxids(
        self, accessions: Dict[str, Set[str]]
    ) -> Tuple[Dict[str, Set[int]], pd.DataFrame]:
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
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--taxa-files",
        help="Path to taxonomy.tsv files for GTDB, either ar53_taxonomy or bac120_taxonomy or both.",
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "--faa-tarball",
        help="Path to tarball containing GTDB ref genome animo acid data sequences",
        required=True,
    )
    parser.add_argument(
        "--dbdir", help="Path to output GTDB database directory", required=True
    )
    parser.add_argument(
        "--tmpdir",
        help="Path to temporary directory",
        required=False,
    )
    parser.add_argument(
        "--keep-temp",
        help="Keep temporary files.",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    make_taxdump_files(args.taxa_files, args.dbdir)

    make_gtdb_db(
        args.faa_tarball,
        os.path.join(args.dbdir, "names.dmp"),
        os.path.join(args.dbdir, "nodes.dmp"),
        args.dbdir,
        args.tmpdir,
        args.keep_temp,
    )


if __name__ == "__main__":
    # python -m autometa.taxonomy.gtdb --dbdir </path/to/gtdb>
    main()
