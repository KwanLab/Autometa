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
params.cov_stddev_limit = 25.0
params.gc_stddev_limit = 5.0
params.interim = "</path/to/store/user/interimediate/results>"
params.processed = "</path/to/store/user/final/results>"


process BINNING {
  tag "Performing Autometa binning"
  container = 'jason-c-kwan/autometa:dev'
  publishDir params.processed, pattern: "${coverage.simpleName}.${params.kingdom}.*.tsv"

  input:
    path kmers
    path coverage
    path gc_content
    path markers
    path taxonomy

  output:
    path "${coverage.simpleName}.${params.kingdom}.binning.tsv", emit: binning
    path "${coverage.simpleName}.${params.kingdom}.kmers.embedded.tsv", emit: embedded_kmers

  """
  autometa-binning \
    --kmers $kmers \
    --coverage $coverage \
    --gc-content $gc_content \
    --markers $markers \
    --output ${coverage.simpleName}.${params.kingdom}.binning.tsv \
    --embedded-kmers ${coverage.simpleName}.${params.kingdom}.kmers.embedded.tsv \
    --embedding-method ${params.kmer_embed_method} \
    --clustering-method ${params.clustering_method} \
    --completeness ${params.completeness} \
    --purity ${params.purity} \
    --cov-stddev-limit ${params.cov_stddev_limit} \
    --gc-stddev-limit ${params.gc_stddev_limit} \
    --taxonomy $taxonomy \
    --starting-rank ${params.binning_starting_rank} \
    --domain ${params.kingdom}
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
    --kmers $kmers \
    --coverage $coverage \
    --binning $assignments \
    --markers $markers \
    --output ${coverage.simpleName}.${params.kingdom}.recruitment.tsv
  """
}
