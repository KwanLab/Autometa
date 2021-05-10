#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process SORT_READS {
  tag "Sorting reads to ${sam.simpleName}"
  cpus params.cpus

  input:
    path sam

  output:
    path "${sam.simpleName}.bam"

  """
  samtools view -@${task.cpus} -bS ${sam} \
    | samtools sort -@${task.cpus} -o ${sam.simpleName}.bam
  """
}
