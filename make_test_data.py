#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

Generate all required intermediate objects and files for Autometa unit testing.
Once generated, each is placed in the `self.data` dict with its respective stage.
After all required objects are within `self.data`, the dictionary is written as
a json file to `test_data.json` to be read by pytest-variables during testing.

`self.data` keys and subkeys for accession with `variables` in test_*.py scripts.

1 metagenome
    1.1 assembly
    1.2 orfs
2 kmers
    2.1 counts
    2.2 am_clr_normalized_counts
3 coverage
    3.1 spades_records
    3.2 sam
    3.3 bam
    3.4 bed
    3.5 fwd_reads
    3.6 rev_reads
4 markers
    4.1 scans
    4.2 filtered_markers
    4.3 orfs
5 taxonomy
    5.1 prot_orfs
    5.2 blastp
    5.3 acc2taxid
    5.4 merged
    5.5 nodes
    5.6 names
6 binning
    6.1 kmers_normalized
    6.2 kmers_embedded
    6.3 taxonomy
    6.4 coverage
    6.5 gc_content
    6.6 markers
7 summary
    7.1 bin_df
"""


import gzip
import json
import os
import typing

import attr
import logging

import pandas as pd
from Bio import SeqIO

from autometa.common import kmers, markers
from autometa.common.external import hmmer, prodigal
from autometa.taxonomy.ncbi import NCBI

import numpy as np

logger = logging.getLogger(__name__)

logging.basicConfig(
    format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
    level=logging.DEBUG,
)


def opener(fpath, mode="r"):
    if fpath.endswith(".gz"):
        if "r" in mode:
            mode += "t"
        fh = gzip.open(fpath, mode)
    else:
        fh = open(fpath, mode)
    return fh


def subset_acc2taxids(blastp_accessions: set, ncbi: NCBI) -> dict:
    acc2taxids = {}
    with opener(ncbi.accession2taxid_fpath) as fh:
        fh.readline()  # skip reading header line
        for line in fh:
            acc_num, acc_ver, taxid, _ = line.split("\t")
            if acc_num in blastp_accessions:
                acc2taxids[acc_num] = taxid
            if acc_ver in blastp_accessions:
                acc2taxids[acc_ver] = taxid
    return acc2taxids


@attr.s(auto_attribs=True)
class TestData:
    metagenome: str
    metagenome_nucl_orfs: str
    metagenome_prot_orfs: str
    coverage_sam: str
    coverage_bed: str
    fwd_reads: str
    rev_read: str
    markers_orfs: str
    markers_scans: str
    markers_filtered: str
    taxonomy_ncbi: str
    taxonmy_blastp: str
    taxonomy_orfs: str
    binning_norm_kmers: str
    binning_embedded_kmers: str
    binning_coverage: str
    binning_gc_content: str
    binning_markers: str
    binning_taxonomy: str
    summary_bin_df: str
    recruitment_binning: str
    data: typing.Dict = attr.ib(default={})
    seed: int = 42

    def prepare_metagenome(self, num_records: int = 4):
        logger.info("Preparing metagenome records test data...")
        records = {}
        for record in SeqIO.parse(self.metagenome, "fasta"):
            records.update({f">{record.id}": str(record.seq)})
            if len(records) >= num_records:
                break

        try:
            prodigal.run(
                assembly=self.metagenome,
                nucls_out=self.metagenome_nucl_orfs,
                prots_out=self.metagenome_prot_orfs,
                force=False,
            )
        except FileExistsError:
            logger.debug("metagenome orfs already exist")
        amino_acid_orfs = {
            f">{orf.id}": str(orf.seq)
            for orf in SeqIO.parse(self.metagenome_prot_orfs, "fasta")
            if f">{orf.id.rsplit('_', 1)[0]}" in records
        }
        self.data["metagenome"] = {"assembly": records, "orfs": amino_acid_orfs}

    def get_kmers(self, num_records: int = 5):
        if num_records < 5:
            raise ValueError(
                f"At least 5 records are required for embedding tests! provided: {num_records}"
            )
        logger.info("Preparing kmer counts test data...")
        # kmer size is 5 (b/c this is the default).
        counts = kmers.count(assembly=self.metagenome, size=5)
        # subset counts to `num_records`
        counts = counts.sample(n=num_records, random_state=42)
        # method is am_clr (b/c this is the default).
        am_clr_normalized_counts = kmers.normalize(df=counts, method="am_clr")

        for df in [counts, am_clr_normalized_counts]:
            df.reset_index(inplace=True)
        self.data["kmers"] = {
            "counts": counts.to_json(),
            "am_clr_normalized_counts": am_clr_normalized_counts.to_json(),
        }

    def get_markers(self):
        logger.info("Preparing orfs for markers annotation")
        try:
            prodigal.run(
                assembly=self.metagenome,
                nucls_out=self.metagenome_nucl_orfs,
                prots_out=self.metagenome_prot_orfs,
                force=False,
            )
        except FileExistsError:
            logger.debug("markers orfs already exist")
        markers_query_orfs = [
            record
            for record in SeqIO.parse(self.metagenome_prot_orfs, "fasta")
            if record.id == "NODE_1505_length_7227_cov_222.087_6"
        ]
        if not os.path.exists(self.markers_orfs):
            SeqIO.write(markers_query_orfs, self.markers_orfs, "fasta")
        markers_query_orfs = {f">{rec.id}": str(rec.seq) for rec in markers_query_orfs}
        logger.info("Annotating ORFs with single-copy markers")
        if not os.path.exists(self.markers_scans) or not os.path.exists(
            self.markers_filtered
        ):
            self.markers_filtered = (
                self.markers_filtered.replace(".gz", "")
                if self.markers_filtered.endswith(".gz")
                else self.markers_filtered
            )
            markers.get(
                kingdom="archaea",
                orfs=self.markers_orfs,
                dbdir=markers.MARKERS_DIR,
                scans=self.markers_scans,
                out=self.markers_filtered,
                seed=self.seed,
            )
        # Retrieve test output hmmscan table
        scans = pd.read_csv(self.markers_scans, sep="\s+", header=None, comment="#")
        filtered_markers = pd.read_csv(self.markers_filtered, sep="\t")
        # The ORFs are necessary for ORF to contig translations
        self.data["markers"] = {
            "scans": scans.to_json(),
            "filtered_markers": filtered_markers.to_json(),
            "orfs": markers_query_orfs,
        }

    def get_taxonomy(self, num_orfs: int = 2):
        logger.info("Making taxonomy test data...")
        # Get diamond blastp output table
        orf_column = 0
        blastp = pd.read_csv(
            self.taxonmy_blastp, sep="\t", index_col=orf_column, header=None
        )
        # Get number of unique ORFs set by `num_orfs`, default is 2.
        orf_hits = set(blastp.index.unique().tolist()[:num_orfs])
        blastp = blastp.loc[orf_hits]
        blastp.reset_index(inplace=True)
        if num_orfs == 2:
            # NODE_38_length_280079_cov_224.186_1 and NODE_38_length_280079_cov_224.186_2
            # together have 400 hits
            assert blastp.shape == (
                400,
                12,
            ), f"shape: {blastp.shape}\ncolumns: {blastp.columns}"

        blastp_query_orfs = {
            f">{record.id}": str(record.seq)
            for record in SeqIO.parse(self.taxonomy_orfs, "fasta")
            if not record.id in orf_hits
        }

        ncbi = NCBI(self.taxonomy_ncbi)
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

    def get_bed_alignments(self, num_contigs=1):
        contig_col = 0
        coverages = pd.read_csv(
            self.coverage_bed, sep="\t", index_col=contig_col, header=None
        )
        # Get number of unique contigs set by `num_contigs`, default is 1.
        coverages = coverages.sample(n=num_contigs, random_state=self.seed)
        coverages.reset_index(inplace=True)
        # Here we are ready to send to json object self.data
        return coverages

    def get_sam_alignments(self, num_contigs=1):
        with opener(self.coverage_sam) as fh:
            lines = ""
            contig = None
            for line in fh:
                if line.startswith("@SQ"):
                    # example line
                    # @SQ	SN:NODE_1503_length_7231_cov_222.076	LN:7231
                    contig = line.split("\t")[1]
                    contig = contig.replace("SN:", "")
                if line.startswith("@HD") or line.startswith("@PG") or contig in line:
                    # @HD and @PG are required for bam construction
                    lines += line
        return lines

    def get_reads(self, read_count=5):
        reads = []
        for file in [self.fwd_reads, self.rev_read]:
            outlines = ""
            count = 0
            with opener(file) as fh:
                for line in fh:
                    if "+" in line:
                        count += 1
                    if count >= read_count:
                        break
                    outlines += line
            reads.append(outlines)
        return reads

    def get_coverage(self):
        logging.info("Making alignment records (sam, and bed files) ...")
        bed = self.get_bed_alignments()
        sam = self.get_sam_alignments()
        logging.info("Getting fwd and rev reads ...")
        fwd_reads, rev_reads = self.get_reads()
        self.data["coverage"] = {
            "bed": bed.to_json(),
            "sam": sam,
            "fwd_reads": fwd_reads,
            "rev_reads": rev_reads,
        }

    def get_binning(self, num_contigs: int = None):
        # Need kmers, coverage, markers, taxonomy
        logger.info("Preparing binning test data")
        annotations = {
            "kmers_normalized": self.binning_norm_kmers,
            "kmers_embedded": self.binning_embedded_kmers,
            "taxonomy": self.binning_taxonomy,
            "coverage": self.binning_coverage,
            "gc_content": self.binning_gc_content,
        }

        markers_df = pd.read_csv(self.binning_markers, sep="\t", index_col="contig")
        contigs = None
        for annotation, fpath in annotations.items():
            df = pd.read_csv(fpath, sep="\t", index_col="contig")
            # We'll grab the first `num_contigs` from the first dataframe (kmers)
            if not contigs and num_contigs:
                # We need to ensure the contigs we pull contain markers...
                contigs = set(
                    df[df.index.isin(markers_df.index)].index.tolist()[:num_contigs]
                )
            if annotation == "taxonomy":
                for column in df.select_dtypes(object).columns:
                    df[column] = df[column].map(lambda taxon: taxon.lower())
            # We need to reset the index from contig to None before json export.
            if contigs:
                jsonified = df.loc[contigs].reset_index().to_json()
            else:
                jsonified = df.reset_index().to_json()
            if "binning" not in self.data:
                self.data["binning"] = {annotation: jsonified}
            else:
                self.data["binning"].update({annotation: jsonified})
        markers_df.reset_index(inplace=True)
        self.data["binning"].update({"markers": markers_df.to_json()})

    def get_summary(self):
        bin_df = pd.read_csv(self.summary_bin_df, sep="\t")
        if "coverage" not in bin_df.columns:
            bin_df["coverage"] = bin_df.contig.map(lambda x: x.split("_cov_")[-1])
        if "GC" not in bin_df.columns:
            bin_df["GC"] = np.random.random_sample(bin_df.contig.nunique())
        if "length" not in bin_df.columns:
            bin_df["length"] = bin_df.contig.map(
                lambda x: x.split("_length_")[-1].split("_cov_")[0]
            )
        self.data["summary"] = {"bin_df": bin_df.to_json()}

    def get_recruitment(self, num_contigs: int = None):
        logger.info("Preparing recruitment binning test data")
        df = pd.read_csv(self.recruitment_binning, sep="\t")
        if num_contigs:
            df = df.sample(n=num_contigs, random_state=self.seed)
        self.data["recruitment"] = {"binning": df.to_json()}

    def to_json(self, out: str):
        logger.info(f"Serializing data to {out}")
        with opener(out, "w") as fh:
            json.dump(obj=self.data, fp=fh)
        logger.info(f"Wrote test data to {out}")


def main():

    outdir = os.path.join("tests", "data")
    metagenome = os.path.join(outdir, "records.fna")
    metagenome_nucl_orfs = os.path.join(outdir, "metagenome_nucl_orfs.fasta")
    metagenome_prot_orfs = os.path.join(outdir, "metagenome_prot_orfs.fasta")
    markers_orfs = os.path.join(outdir, "markers_orfs.faa")
    markers_scans = os.path.join(outdir, "markers_scans.tsv.gz")
    markers_filtered = os.path.join(outdir, "markers_filtered.tsv.gz")
    taxonomy_ncbi = os.path.join("autometa", "databases", "ncbi")
    taxonmy_blastp = os.path.join(outdir, "blastp.tsv.gz")
    taxonomy_orfs = os.path.join(outdir, "taxonomy_orfs.faa")
    binning_norm_kmers = os.path.join(outdir, "binning_kmers.am_clr.tsv.gz")
    binning_embedded_kmers = os.path.join(outdir, "binning_kmers.am_clr.bhsne.tsv.gz")
    binning_coverage = os.path.join(outdir, "binning_coverage.tsv.gz")
    binning_gc_content = os.path.join(outdir, "binning_gc_content.tsv.gz")
    binning_markers = os.path.join(outdir, "binning_markers.tsv.gz")
    binning_taxonomy = os.path.join(outdir, "binning_taxonomy.tsv.gz")
    coverage_sam = os.path.join(outdir, "records.sam.gz")
    coverage_bed = os.path.join(outdir, "records.bed.gz")
    fwd_reads = os.path.join(outdir, "records_1.fastq.gz")
    rev_read = os.path.join(outdir, "records_2.fastq.gz")
    summary_bin_df = os.path.join(outdir, "summary_bin_df.tsv.gz")
    recruitment_binning = os.path.join(outdir, "recruitment_binning.tsv.gz")

    test_data = TestData(
        metagenome=metagenome,
        metagenome_nucl_orfs=metagenome_nucl_orfs,
        metagenome_prot_orfs=metagenome_prot_orfs,
        markers_orfs=markers_orfs,
        markers_scans=markers_scans,
        markers_filtered=markers_filtered,
        taxonomy_ncbi=taxonomy_ncbi,
        taxonmy_blastp=taxonmy_blastp,
        taxonomy_orfs=taxonomy_orfs,
        binning_norm_kmers=binning_norm_kmers,
        binning_embedded_kmers=binning_embedded_kmers,
        binning_coverage=binning_coverage,
        binning_gc_content=binning_gc_content,
        binning_markers=binning_markers,
        binning_taxonomy=binning_taxonomy,
        coverage_sam=coverage_sam,
        coverage_bed=coverage_bed,
        fwd_reads=fwd_reads,
        rev_read=rev_read,
        summary_bin_df=summary_bin_df,
        recruitment_binning=recruitment_binning,
    )

    # TODO: Decrease the size of the test_data.json file...
    test_data.prepare_metagenome()
    test_data.get_kmers()
    # COMBAK: Minimize data structures for coverage test data
    test_data.get_coverage()
    # # COMBAK: Minimize data structures for taxonomy test data
    test_data.get_taxonomy()
    test_data.get_markers()
    test_data.get_binning(num_contigs=0)  # all contigs
    test_data.get_recruitment()
    test_data.get_summary()

    out = os.path.join(outdir, "test_data.json")
    test_data.to_json(out=out)
    return


if __name__ == "__main__":
    main()
