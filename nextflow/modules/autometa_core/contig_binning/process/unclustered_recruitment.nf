#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process UNCLUSTERED_RECRUITMENT {
  tag "Performing Autometa unclustered recruitment"
  publishDir params.outdir, pattern: "${coverage.simpleName}.${params.kingdom}.recruitment.tsv.gz"

  input:
    path kmers
    path coverage
    path binning
    path markers
    path taxonomy

  output:
    path "${coverage.simpleName}.${params.kingdom}.recruitment.tsv.gz", emit: binning
    path "${coverage.simpleName}.${params.kingdom}.recruitment.main.tsv.gz", emit: main

  """
  autometa-unclustered-recruitment \
    --classifier ${params.classification_method} \
    --kmer-dimensions ${params.classification_kmer_pca_dimensions} \
    --seed 42 \
    --taxonomy $taxonomy \
    --kmers $kmers \
    --coverage $coverage \
    --binning $binning \
    --markers $markers \
    --output-binning ${coverage.simpleName}.${params.kingdom}.recruitment.tsv.gz \
    --output-main ${coverage.simpleName}.${params.kingdom}.recruitment.main.tsv.gz
  """
}
