#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ORFS {
  tag "Calling orfs for ${metagenome.simpleName}"
  // Hardcoding cpus here b/c prodigal is limited to only using single core
  cpus 1
  publishDir params.interim_dir, pattern: "${metagenome.simpleName}.orfs.f*"

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.orfs.fna", emit: nucls
    path "${metagenome.simpleName}.orfs.faa", emit: prots

  """
  prodigal \
    -i $metagenome \
    -d ${metagenome.simpleName}.orfs.fna \
    -a ${metagenome.simpleName}.orfs.faa \
    -p meta \
    -q \
    -m
  """
}