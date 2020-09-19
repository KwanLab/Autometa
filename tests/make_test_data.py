import json

import pandas as pd

from Bio import SeqIO


from autometa.common import kmers, markers
from autometa.common.external import prodigal


# TODO: Decrease the size of the test_data.json file by only using
# a few contigs rather than currently the 50 in records.fna

# TODO: Will need to replace hardcoded records.fna with
# from autometa.datasets import load_test_data
# test_data_fpath = load_test_data()
# counts = kmers.count(assembly=test_data_fpath, size=5) ... etc...


counts = kmers.count(assembly="records.fna", size=5)
norm_df = kmers.normalize(df=counts)

records = {
    f">{record.id}": str(record.seq) for record in SeqIO.parse("records.fna", "fasta")
}
try:
    prodigal.run(
        assembly="records.fna", nucls_out="orfs.fna", prots_out="orfs.faa", force=True
    )
except FileExistsError:
    pass


orfs = {f">{record.id}": str(record.seq) for record in SeqIO.parse("orfs.faa", "fasta")}

bact_markers = markers.get(
    orfs="orfs.faa", kingdom="bacteria", dbdir=markers.MARKERS_DIR
)
arch_markers = markers.get(
    orfs="orfs.faa", kingdom="archaea", dbdir=markers.MARKERS_DIR
)


# counts.reset_index(inplace=True)
# norm_df.reset_index(inplace=True)
# bact_markers.reset_index(inplace=True)
# arch_markers.reset_index(inplace=True)
for df in [counts, norm_df, bact_markers, arch_markers]:
    df.reset_index(inplace=True)


d = {
    "kmers": {
        "counts": counts.to_json(),
        "norm_df": norm_df.to_json(),
        "small_metagenome": records,
    },
    "coverage": {
        "spades_records": "records.fna",
        "bam": "alignments.bam",
        "sam": "alignments.sam",
    },
    "markers": {
        "orfs": orfs,
        "bacteria": bact_markers.to_json(),
        "archaea": arch_markers.to_json(),
    },
}


with open("test_data.json", "w") as fh:
    json.dump(obj=d, fp=fh)
