#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process ANALYZE_KMERS {
  tag "counting kmers for ${metagenome.simpleName}"
  label "python_cpus"
  publishDir params.interim_dir, pattern: "*.kmers.*"

  input:
    path metagenome

  output:
    path "*.kmers.tsv", emit: counts
    path "*.kmers.normalized.tsv", emit: normalized
    path "*.kmers.embedded.tsv", emit: embedded

  """
  autometa-kmers \
    --fasta $metagenome \
    --kmers ${metagenome.simpleName}.kmers.tsv \
    --size ${params.kmer_size} \
    --norm-output ${metagenome.simpleName}.kmers.normalized.tsv \
    --norm-method ${params.norm_method} \
    --pca-dimensions ${params.pca_dimensions} \
    --embedding-output ${metagenome.simpleName}.kmers.embedded.tsv \
    --embedding-method ${params.embedding_method} \
    --embedding-dimensions ${params.embedding_dimensions} \
    --cpus ${task.cpus} \
    --seed 42
  """
}