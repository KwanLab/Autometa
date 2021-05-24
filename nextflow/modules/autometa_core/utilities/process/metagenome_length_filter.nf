#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process LENGTH_FILTER {
  tag "filtering metagenome ${metagenome.simpleName}"
  publishDir params.interim_dir, pattern: "${metagenome.simpleName}.*"

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.filtered.fna", emit: fasta
    path "${metagenome.simpleName}.stats.tsv", emit: stats
    path "${metagenome.simpleName}.gc_content.tsv", emit: gc_content

  script:
  if ( "${metagenome.Name}" == "${metagenome.simpleName}.filtered.fna" )
  """
  autometa-length-filter \
    --assembly $metagenome \
    --cutoff ${params.length_cutoff} \
    --output-fasta "${metagenome.simpleName}.autometa_filtered.fna" \
    --output-stats ${metagenome.simpleName}.stats.tsv \
    --output-gc-content ${metagenome.simpleName}.gc_content.tsv
  """
  else 
  """
  autometa-length-filter \
    --assembly $metagenome \
    --cutoff ${params.length_cutoff} \
    --output-fasta "${metagenome.simpleName}.filtered.fna" \
    --output-stats ${metagenome.simpleName}.stats.tsv \
    --output-gc-content ${metagenome.simpleName}.gc_content.tsv
  """
}
