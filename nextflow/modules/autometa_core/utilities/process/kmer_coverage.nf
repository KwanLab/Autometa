#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SPADES_KMER_COVERAGE {
  tag "Calculating k-mer coverage for ${metagenome.simpleName}"
  
  publishDir params.interim_dir, pattern: "${metagenome.simpleName}.coverages.tsv"

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.coverages.tsv"

  """
  autometa-coverage \
    --assembly $metagenome \
    --from-spades \
    --out ${metagenome.simpleName}.coverages.tsv
  """
}
