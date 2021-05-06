#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process LCA {
  label 'process_medium'
  label 'process_long'
  
  tag "Assigning LCA to ${blast.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:rw"
  publishDir params.interim_dir, pattern: "${blast.simpleName}.lca.tsv"

  input:
    path blast

  output:
    path "${blast.simpleName}.lca.tsv"

  """
  autometa-taxonomy-lca --blast ${blast} --dbdir /ncbi --output ${blast.simpleName}.lca.tsv
  """
}