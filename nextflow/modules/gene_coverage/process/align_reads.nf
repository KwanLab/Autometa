#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process ALIGN_READS {
  tag "Aligning reads to ${metagenome.simpleName}"
  label "python_cpus"

  input:
    path metagenome
    path fwd_reads
    path rev_reads
    path se_reads

  output:
    path "${metagenome.simpleName}.sam"

  """
  bowtie2-build \
    --threads ${task.cpus}
    ${metagenome} \
    ${metagenome.simpleName}.db

  bowtie2 \
    -x ${metagenome.simpleName}.db \
    -q \
    --phred33 \
    --very-sensitive \
    --no-unal \
    -p ${task.cpus} \
    -S ${metagenome.simpleName}.sam \
    -1 $fwd_reads \
    -2 $rev_reads \
    -U $se_reads
  """
}
