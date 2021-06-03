#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process KMER_COVERAGE {
  tag "Retrieving k-mer coverage from ${metagenome.simpleName} headers"
  label "python_cpus"
  
  publishDir params.interim_dir, pattern: "${metagenome.simpleName}.coverages.tsv"

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.coverages.tsv"

  """
  autometa-coverage \
    --assembly $metagenome \
    --cpus ${task.cpus} \
    --from-spades \
    --out ${metagenome.simpleName}.coverages.tsv
  """
}
