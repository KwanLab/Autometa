#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MAJORITY_VOTE {
  label 'process_medium'
  
  tag "Performing taxon majority vote on ${lca.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:rw"
  publishDir params.interim_dir, pattern: "${lca.simpleName}.votes.tsv"

  input:
    path lca

  output:
    path "${lca.simpleName}.votes.tsv"

  """
  autometa-taxonomy-majority-vote --lca ${lca} --output ${lca.simpleName}.votes.tsv --dbdir /ncbi
  """
}