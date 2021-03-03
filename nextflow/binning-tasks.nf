#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.kingdom = "bacteria"
params.completeness = 20.0
params.purity = 90.0
params.interim = "</path/to/store/user/interimediate/results>"
params.processed = "</path/to/store/user/final/results>"


process BINNING {
  tag "Performing Autometa binning"
  // container = 'placeholder for autometa image'
  publishDir params.processed, pattern: "${coverage.simpleName}.${params.kingdom}.*.tsv", mode:'copy'

  input:
    path kmers
    path coverage
    path markers
    path taxonomy

  output:
    path "${coverage.simpleName}.${params.kingdom}.binning.tsv", emit: binning
    path "${coverage.simpleName}.${params.kingdom}.kmers.embedded.tsv", emit: embedded_kmers

  """
  autometa-binning \
    --clustering-method dbscan \
    --completeness ${params.completeness} \
    --purity ${params.purity} \
    --taxonomy $taxonomy \
    --starting-rank superkingdom \
    --domain ${params.kingdom} \
    --embedding-method bhsne \
    --embedded-kmers ${coverage.simpleName}.${params.kingdom}.kmers.embedded.tsv \
    $kmers \
    $coverage \
    $markers \
    ${coverage.simpleName}.${params.kingdom}.binning.tsv
  """
}

process UNCLUSTERED_RECRUITMENT {
  tag "Performing Autometa unclustered recruitment"
  // container = 'placeholder for autometa image'
  publishDir params.processed, pattern: "${coverage.simpleName}.${params.kingdom}.recruitment.tsv", mode:'copy'

  input:
    path kmers
    path coverage
    path assignments
    path markers
    path taxonomy

  output:
    path "${coverage.simpleName}.${params.kingdom}.recruitment.tsv", emit: binning

  """
  autometa-unclustered-recruitment \
    --classifier decision_tree \
    --kmer-dimensions 50 \
    --seed 42 \
    --taxonomy $taxonomy \
    $kmers \
    $coverage \
    $assignments \
    $markers \
    ${coverage.simpleName}.${params.kingdom}.recruitment.tsv
  """
}
