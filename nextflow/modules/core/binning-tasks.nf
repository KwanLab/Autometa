#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Data inputs
params.interim_dir = "</path/to/store/user/interimediate/results>"
params.outdir = "</path/to/store/user/final/results>"
// Binning parameters
params.kingdom = "bacteria"
params.clustering_method = "dbscan" // choices: "dbscan", "hdbscan"
params.binning_starting_rank = "superkingdom" // choices: "superkingdom", "phylum", "class", "order", "family", "genus", "species"
params.completeness = 20.0
params.purity = 95.0
params.cov_stddev_limit = 25.0
params.gc_stddev_limit = 5.0
// Unclustered recruitment parameters
params.classification_method = "decision_tree" // choices: "decision_tree", "random_forest"
params.classification_kmer_pca_dimensions = 50
// Summary parameters
params.ncbi_database = "Path to user ncbi databases directory"

process BINNING {
  tag "Performing Autometa binning"
  container = 'chaseauto:latest'
  publishDir params.outdir, pattern: "${coverage.simpleName}.${params.kingdom}.*.tsv.gz"

  input:
    path kmers
    path coverage
    path gc_content
    path markers
    path taxonomy

  output:
    path "${coverage.simpleName}.${params.kingdom}.binning.tsv.gz", emit: binning
    path "${coverage.simpleName}.${params.kingdom}.main.tsv.gz", emit: main

  """
  autometa-binning \
    --kmers $kmers \
    --coverages $coverage \
    --gc-content $gc_content \
    --markers $markers \
    --output-binning ${coverage.simpleName}.${params.kingdom}.binning.tsv.gz \
    --output-main ${coverage.simpleName}.${params.kingdom}.main.tsv.gz \
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
  container = 'chaseauto:latest'
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

process BINNING_SUMMARY {
  tag "Binning summary for ${binning_main.simpleName}"
  container = 'chaseauto:latest'
  containerOptions = "-v ${params.single_db_dir}:/ncbi:ro"

  input:
    path binning_main
    path markers
    path metagenome
    val binning_column

  output:
    path 'metabin_stats.tsv', emit: stats
    path 'metabin_taxonomy.tsv', emit: taxonomies
    path 'metabins', emit: metabins

  script:
  """
  autometa-binning-summary \
    --ncbi /ncbi \
    --binning-main $binning_main \
    --markers $markers \
    --metagenome $metagenome \
    --binning-column $binning_column \
    --output-stats metabin_stats.tsv \
    --output-taxonomy metabin_taxonomy.tsv \
    --output-metabins metabins
  """
}
