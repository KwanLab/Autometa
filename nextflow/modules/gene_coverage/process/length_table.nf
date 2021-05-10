#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LENGTH_TABLE {
  tag "length table for ${metagenome.simpleName}"
  cpus params.cpus

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.lengths.tsv"

  """
  #!/usr/bin/env python

  from Bio import SeqIO
  import pandas as pd

  seqs = {record.id: len(record) for record in SeqIO.parse($metagenome, "fasta")}
  lengths = pd.Series(seqs, name="length")
  lengths.index.name = "contig"
  lengths.to_csv(${metagenome.simpleName}.lengths.tsv, sep="\t", index=True, header=True)
  """
}