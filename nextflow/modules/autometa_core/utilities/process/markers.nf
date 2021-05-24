#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MARKERS {
  tag "Finding markers for ${orfs.simpleName}"
  label "process_low"
  
  // copying orfs via stageInMode is required to run hmmscan (does not handle symlinks)
  stageInMode 'copy'
  publishDir params.interim_dir, pattern: "${orfs.simpleName}.markers.tsv"
  publishDir params.interim_dir, pattern: "${orfs.simpleName}.hmmscan.tsv"

  input:
    path orfs

  output:
    path "${orfs.simpleName}.markers.tsv"

  """
  autometa-markers \
    --orfs $orfs \
    --hmmscan ${orfs.simpleName}.hmmscan.tsv \
    --out ${orfs.simpleName}.markers.tsv \
    --kingdom ${params.kingdom} \
    --parallel \
    --cpus ${task.cpus} \
    --seed 42
  """
}
