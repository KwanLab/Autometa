#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.kingdom = "bacteria"
params.classification_kmer_pca_dimensions = 50
params.kmer_embed_method = "bhsne" // choices: "bhsne", "sksne", "umap"
params.clustering_method = "dbscan" // choices: "dbscan", "hdbscan"
params.binning_starting_rank = "superkingdom" // choices: "superkingdom", "phylum", "class", "order", "family", "genus", "species"
params.classification_method = "decision_tree" // choices: "decision_tree", "random_forest"
params.completeness = 20.0
params.purity = 90.0
params.interim = "</path/to/store/user/interimediate/results>"
params.processed = "</path/to/store/user/final/results>"


process BINNING {
  tag "Performing Autometa binning"
  container = 'jason-c-kwan/autometa:dev'
  publishDir params.processed, pattern: "${coverage.simpleName}.${params.kingdom}.*.tsv"

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
    --clustering-method ${params.clustering_method} \
    --completeness ${params.completeness} \
    --purity ${params.purity} \
    --taxonomy $taxonomy \
    --starting-rank ${params.binning_starting_rank} \
    --domain ${params.kingdom} \
    --embedding-method ${params.kmer_embed_method} \
    --embedded-kmers ${coverage.simpleName}.${params.kingdom}.kmers.embedded.tsv \
    $kmers \
    $coverage \
    $markers \
    ${coverage.simpleName}.${params.kingdom}.binning.tsv
  """
}

process UNCLUSTERED_RECRUITMENT {
  tag "Performing Autometa unclustered recruitment"
  container = 'jason-c-kwan/autometa:dev'
  publishDir params.processed, pattern: "${coverage.simpleName}.${params.kingdom}.recruitment.tsv"

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
    --classifier ${params.classification_method} \
    --kmer-dimensions ${params.classification_kmer_pca_dimensions} \
    --seed 42 \
    --taxonomy $taxonomy \
    $kmers \
    $coverage \
    $assignments \
    $markers \
    ${coverage.simpleName}.${params.kingdom}.recruitment.tsv
  """
}
