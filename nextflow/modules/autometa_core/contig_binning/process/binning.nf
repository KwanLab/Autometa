#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BINNING {
  tag "Performing Autometa binning"
  publishDir params.outdir, pattern: "${coverage.simpleName}.${params.kingdom}.*.tsv.gz"

  input:
    path kmers
    path coverage
    path gc_content
    path markers
    val taxonomy

  output:
    path "${coverage.simpleName}.${params.kingdom}.binning.tsv.gz", emit: binning
    path "${coverage.simpleName}.${params.kingdom}.main.tsv.gz", emit: main

  script:
  if (taxonomy == false) 
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
        --starting-rank ${params.binning_starting_rank} \
        --domain ${params.kingdom}
      """
  else
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
