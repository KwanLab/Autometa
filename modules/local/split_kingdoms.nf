#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SPLIT_KINGDOMS {
  label 'process_medium'
  
  tag "Splitting votes into kingdoms for ${assembly.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:rw"
  publishDir params.interim_dir, pattern: "${assembly.simpleName}.taxonomy.tsv"
  publishDir params.interim_dir, pattern: '*.{bacteria,archaea}.fna'

  input:
    path votes
    path assembly

  output:
    path "${assembly.simpleName}.taxonomy.tsv", emit: taxonomy
    path "${assembly.simpleName}.bacteria.fna", emit: bacteria
    path "${assembly.simpleName}.archaea.fna", emit: archaea, optional:true

  """
  autometa-taxonomy \
    --votes ${votes} \
    --output . \
    --prefix ${assembly.simpleName} \
    --split-rank-and-write superkingdom \
    --assembly ${assembly} \
    --ncbi /ncbi
  """
}
