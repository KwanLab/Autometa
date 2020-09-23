import json

import pandas as pd

from Bio import SeqIO


from autometa.common import kmers, markers
from autometa.common.external import prodigal
from autometa.taxonomy.ncbi import NCBI
import gzip
import attr
from pathlib import Path


def subset_acc2taxids(blastp: str, ncbi: NCBI) -> dict:
    accessions = set()
    add_accessions = accessions.add
    with open("tests/data/metagenome.filtered.orfs.dmnd.blastp") as fh:
        for line in fh:
            subject_accession = line.strip().split("\t")[1]
            add_accessions(subject_accession)
    acc2taxids = {}
    fh = (
        gzip.open(ncbi.accession2taxid_fpath, "rt")
        if ncbi.accession2taxid_fpath.endswith(".gz")
        else open(ncbi.accession2taxid_fpath)
    )
    for line in fh:
        acc_num, acc_ver, taxid, _ = line.split("\t")
        if acc_num in accessions or acc_ver in accessions:
            acc2taxids.update({acc_num: {"taxid": taxid, "acc_ver": acc_ver}})
    return acc2taxids


@attr.s(kw_only=True)
class TestData:
    records_fpath = attr.ib(validator=attr.validators.instance_of(str))
    nucls_out = attr.ib(validator=attr.validators.instance_of(str))
    prots_out = attr.ib(validator=attr.validators.instance_of(str))
    ncbi_dirpath = attr.ib(validator=attr.validators.instance_of(str))
    blastp_fpath = attr.ib(validator=attr.validators.instance_of(str))
    blastp_orfs_fpath = attr.ib(validator=attr.validators.instance_of(str))
    data = attr.ib(factory=dict)

    def get_kmers(self):
        print("Preparing kmer counts test data...")
        counts = kmers.count(assembly=self.records_fpath, size=5)
        norm_df = kmers.normalize(df=counts)
        print("Preparing metagenome records test data...")
        records = {
            f">{record.id}": str(record.seq)
            for record in SeqIO.parse(self.records_fpath, "fasta")
        }
        for df in [counts, norm_df]:
            df.reset_index(inplace=True)
        self.data["kmers"] = {
            "counts": counts.to_json(),
            "norm_df": norm_df.to_json(),
            "small_metagenome": records,
        }

    def get_markers(self):
        print("Preparing orf calling test data...")
        try:
            prodigal.run(
                assembly=self.records_fpath,
                nucls_out=self.nucls_out,
                prots_out=self.prots_out,
                force=True,
            )
        except FileExistsError:
            pass

        orfs = {
            f">{record.id}": str(record.seq)
            for record in SeqIO.parse(self.prots_out, "fasta")
        }
        print("Preparing markers test data...")
        bact_markers = markers.get(
            orfs=self.prots_out, kingdom="bacteria", dbdir=markers.MARKERS_DIR
        )
        arch_markers = markers.get(
            orfs=self.prots_out, kingdom="archaea", dbdir=markers.MARKERS_DIR
        )

        for df in [bact_markers, arch_markers]:
            df.reset_index(inplace=True)

        self.data["markers"] = {
            "orfs": orfs,
            "bacteria": bact_markers.to_json(),
            "archaea": arch_markers.to_json(),
        }

    def get_taxonomy(self, num_blastp_hits=2):
        print("Making taxonomy test data...")
        ncbi = NCBI(self.ncbi_dirpath)
        # Get nodes.dmp, names.dmp and merged.dmp data structures.
        nodes = ncbi.nodes
        merged = ncbi.merged
        names = ncbi.names
        # Get diamond blastp output table
        blastp = pd.read_csv(self.blastp_fpath, sep="\t", header=None)
        blastp_hits = set(blastp[0].unique().tolist()[:num_blastp_hits])
        blastp = blastp.set_index(0).loc[blastp_hits].reset_index()
        # COMBAK: Need to ensure the blastp ORFs get parsed proper for LCA
        blastp_orfs = {
            f">{record.id}": str(record.seq)
            for record in SeqIO.parse(self.blastp_orfs_fpath, "fasta")
            if not record.id in blastp_hits
        }

        # Get prot.accession2taxid datastructure and subset by taxids encountered in blastp output.
        acc2taxids = subset_acc2taxids(self.blastp_fpath, ncbi)

        self.data["taxonomy"] = {
            "merged": merged,
            "nodes": nodes,
            "names": names,
            "acc2taxid": acc2taxids,
            "prot_orfs": blastp_orfs,
            "blastp": blastp.to_json(),
        }

    def get_coverage(self):
        self.data["coverage"] = {
            "spades_records": self.records_fpath,
            "bam": "alignments.bam",
            "sam": "alignments.sam",
        }

    def to_json(self, out: str):
        print(f"Serializing data to {out}")
        with open(out, "w") as fh:
            json.dump(obj=self.data, fp=fh)
        print(f"Wrote test data to {out}")


def main():

    # outdir = Path("tests/data")
    records_fpath = "tests/data/records.fna"
    # records_fpath = outdir / "records.fna"
    nucls_out = "tests/data/orfs.fna"
    # nucls_out = outdir / "orfs.fna"
    prots_out = "tests/data/orfs.faa"
    # prots_out = outdir / "orfs.faa"
    ncbi_dirpath = "autometa/databases/ncbi/"
    # ncbi_dirpath = "autometa/databases/ncbi/"
    blastp_fpath = "tests/data/metagenome.filtered.orfs.dmnd.blastp"
    # blastp_fpath = outdir / "metagenome.filtered.orfs.dmnd.blastp"
    blastp_orfs_fpath = "tests/data/metagenome.filtered.orfs.faa"
    # blastp_orfs_fpath = outdir / "metagenome.filtered.orfs.faa"

    test_data = TestData(
        records_fpath=records_fpath,
        nucls_out=nucls_out,
        prots_out=prots_out,
        ncbi_dirpath=ncbi_dirpath,
        blastp_fpath=blastp_fpath,
        blastp_orfs_fpath=blastp_orfs_fpath,
    )

    # TODO: Decrease the size of the test_data.json file by only using
    # as few contigs as possible for each test case...
    test_data.get_kmers()
    test_data.get_coverage()
    test_data.get_taxonomy()
    test_data.get_markers()

    out = "tests/data/test_data.json"
    # out = outdir / "test_data.json"
    test_data.to_json(out=out)
    return


if __name__ == "__main__":
    main()
