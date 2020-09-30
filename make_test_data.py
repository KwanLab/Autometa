#!/usr/bin/env python


import gzip
import json

import attr
import pandas as pd
from Bio import SeqIO

from autometa.common import kmers, markers
from autometa.common.external import hmmer, prodigal
from autometa.taxonomy.ncbi import NCBI
import logging
import os

logger = logging.getLogger(__name__)

logging.basicConfig(
    format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)


def subset_acc2taxids(blastp_accessions: set, ncbi: NCBI) -> dict:
    acc2taxids = {}
    fh = (
        gzip.open(ncbi.accession2taxid_fpath, "rt")
        if ncbi.accession2taxid_fpath.endswith(".gz")
        else open(ncbi.accession2taxid_fpath)
    )
    fh.readline()  # skip reading header line
    for line in fh:
        acc_num, acc_ver, taxid, _ = line.split("\t")
        if acc_num in blastp_accessions:
            acc2taxids[acc_num] = taxid
        if acc_ver in blastp_accessions:
            acc2taxids[acc_ver] = taxid
    return acc2taxids


@attr.s(kw_only=True)
class TestData:
    records_fpath = attr.ib(validator=attr.validators.instance_of(str))
    nucls_out = attr.ib(validator=attr.validators.instance_of(str))
    prots_out = attr.ib(validator=attr.validators.instance_of(str))
    markers_query_orfs_fpath = attr.ib(validator=attr.validators.instance_of(str))
    scans_fpath = attr.ib(validator=attr.validators.instance_of(str))
    markers_fpath = attr.ib(validator=attr.validators.instance_of(str))
    ncbi_dirpath = attr.ib(validator=attr.validators.instance_of(str))
    blastp_fpath = attr.ib(validator=attr.validators.instance_of(str))
    blastp_query_orfs_fpath = attr.ib(validator=attr.validators.instance_of(str))
    data = attr.ib(factory=dict)
    seed = attr.ib(default=42)

    def get_kmers(self, num_records=4):
        if num_records < 4:
            raise ValueError(
                f"At least 4 records are required for embedding tests! provided: {num_records}"
            )
        logging.info("Preparing kmer counts test data...")
        # kmer size is 5 (b/c this is the default).
        counts = kmers.count(assembly=self.records_fpath, size=5)
        # subset counts to first rows via `num_records`
        counts = counts.iloc[:num_records]
        # method is am_clr (b/c this is the default).
        am_clr_normalized_counts = kmers.normalize(df=counts, method="am_clr")
        logging.info("Preparing metagenome records test data...")
        records = {}
        for record in SeqIO.parse(self.records_fpath, "fasta"):
            records.update({f">{record.id}": str(record.seq)})
            if len(records) >= num_records:
                break

        for df in [counts, am_clr_normalized_counts]:
            df.reset_index(inplace=True)
        self.data["kmers"] = {
            "counts": counts.to_json(),
            "am_clr_normalized_counts": am_clr_normalized_counts.to_json(),
        }
        self.data["metagenome"] = {"assembly": records}

    def get_markers(self):
        logging.info("Preparing orfs for markers annotation")
        try:
            prodigal.run(
                assembly=self.records_fpath,
                nucls_out=self.nucls_out,
                prots_out=self.prots_out,
                force=True,
            )
        except FileExistsError:
            logging.debug("markers orfs already exist")
        markers_query_orfs = [
            record
            for record in SeqIO.parse(self.prots_out, "fasta")
            if record.id == "NODE_1505_length_7227_cov_222.087_6"
        ]
        if not os.path.exists(self.markers_query_orfs_fpath):
            SeqIO.write(markers_query_orfs, self.markers_query_orfs_fpath, "fasta")
        markers_query_orfs = {f">{rec.id}": str(rec.seq) for rec in markers_query_orfs}
        logging.info("Annotating ORFs with single-copy markers")
        if not os.path.exists(self.scans_fpath) or not os.path.exists(
            self.markers_fpath
        ):
            markers.get(
                kingdom="archaea",
                orfs=self.markers_query_orfs_fpath,
                dbdir=markers.MARKERS_DIR,
                scans=self.scans_fpath,
                out=self.markers_fpath,
                seed=self.seed,
            )
        # Retrieve test output hmmscan table
        scans = pd.read_csv(self.scans_fpath, sep="\s+", header=None, comment="#")
        filtered_markers = pd.read_csv(self.markers_fpath, sep="\t")
        # The ORFs are necessary for ORF to contig translations
        self.data["markers"] = {
            "scans": scans.to_json(),
            "filtered_markers": filtered_markers.to_json(),
            "orfs": markers_query_orfs,
        }

    def get_taxonomy(self, num_orfs=2):
        logging.info("Making taxonomy test data...")
        # Get diamond blastp output table
        blastp = pd.read_csv(self.blastp_fpath, sep="\t", header=None)
        orf_column = 0
        # Get number of unique ORFs set by `num_orfs`, default is 2.
        orf_hits = set(blastp[orf_column].unique().tolist()[:num_orfs])
        blastp = blastp.set_index(orf_column).loc[orf_hits].reset_index()
        if num_orfs == 2:
            # NODE_38_length_280079_cov_224.186_1 and NODE_38_length_280079_cov_224.186_2
            # together have 400 hits
            assert blastp.shape == (
                400,
                12,
            ), f"shape: {blastp.shape}\ncolumns: {blastp.columns}"

        blastp_query_orfs = {
            f">{record.id}": str(record.seq)
            for record in SeqIO.parse(self.blastp_query_orfs_fpath, "fasta")
            if not record.id in orf_hits
        }

        ncbi = NCBI(self.ncbi_dirpath)
        # Get prot.accession2taxid datastructure and subset by taxids encountered in blastp output.
        sacc_column = 1
        blastp_accessions = set(blastp[sacc_column].unique().tolist())
        acc2taxids = subset_acc2taxids(blastp_accessions, ncbi)
        accessions = {k for k in acc2taxids.keys()}
        blastp = blastp.set_index(sacc_column).loc[accessions].reset_index()
        blastp = blastp.set_index(orf_column).reset_index()
        assert blastp.shape[0] == len(
            acc2taxids
        ), f"blastp shape: {blastp.shape}\tnum. acc2taxids: {len(acc2taxids)}"
        # Get nodes.dmp, names.dmp and merged.dmp data structures.
        nodes = ncbi.nodes
        names = ncbi.names
        # Merged are only necessary if taxids have been deprecated or suppressed
        blastp_taxids = acc2taxids.values()
        merged = {old: new for old, new in ncbi.merged.items() if old in blastp_taxids}

        self.data["taxonomy"] = {
            "prot_orfs": blastp_query_orfs,
            "blastp": blastp.to_json(),
            "acc2taxid": acc2taxids,
            "merged": merged,
            "nodes": nodes,
            "names": names,
        }

    def get_coverage(self):
        self.data["coverage"] = {
            "spades_records": self.records_fpath,
            "bam": "alignments.bam",
            "sam": "alignments.sam",
        }

    def to_json(self, out: str):
        logging.info(f"Serializing data to {out}")
        with open(out, "w") as fh:
            json.dump(obj=self.data, fp=fh)
        logging.info(f"Wrote test data to {out}")


def main():

    outdir = os.path.join("tests", "data")
    records_fpath = os.path.join(outdir, "records.fna")
    nucls_out = os.path.join(outdir, "orfs.fna")
    prots_out = os.path.join(outdir, "orfs.faa")
    markers_query_orfs_fpath = os.path.join(outdir, "markers_query_orfs.faa")
    scans_fpath = os.path.join(outdir, "archaea.hmmscan.tsv")
    markers_fpath = os.path.join(outdir, "archaea.markers.tsv")
    ncbi_dirpath = os.path.join("autometa", "databases", "ncbi")
    blastp_fpath = os.path.join(outdir, "metagenome.filtered.orfs.dmnd.blastp")
    blastp_query_orfs_fpath = os.path.join(outdir, "metagenome.filtered.orfs.faa")

    test_data = TestData(
        records_fpath=records_fpath,
        nucls_out=nucls_out,
        prots_out=prots_out,
        markers_query_orfs_fpath=markers_query_orfs_fpath,
        scans_fpath=scans_fpath,
        markers_fpath=markers_fpath,
        ncbi_dirpath=ncbi_dirpath,
        blastp_fpath=blastp_fpath,
        blastp_query_orfs_fpath=blastp_query_orfs_fpath,
    )

    # TODO: Decrease the size of the test_data.json file...
    test_data.get_kmers()
    # COMBAK: Minimize data structures for coverage test data
    test_data.get_coverage()
    # # COMBAK: Minimize data structures for taxonomy test data
    test_data.get_taxonomy()
    test_data.get_markers()

    out = os.path.join(outdir, "test_data.json")
    test_data.to_json(out=out)
    return


if __name__ == "__main__":
    main()
