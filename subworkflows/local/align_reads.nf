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

  bowtie2 \\
      -x ${metagenome.simpleName}.db \\
      -q \\
      --phred33 \\
      --very-sensitive \\
      --no-unal \\
      -p ${task.cpus} \\
      -S ${metagenome.simpleName}.sam \\
      -1 $fwd_reads \\
      -2 $rev_reads \\
      -U $se_read
  """
}

params.bedtools_genomecov_options      = [:]

include { BOWTIE2_ALIGN } from './../../modules/nf-core/software/bedtools/genomecov.nf' addParams( options: params.bedtools_genomecov_options )
include { BEDTOOLS_GENOMECOV } from './../../modules/nf-core/software/bedtools/genomecov.nf' addParams( options: params.bedtools_genomecov_options )

workflow ALIGN_READS {
    take:
    metagenome        // channel: [ val(meta), path(bam) ]
    reads             // channel: [ val(meta), path(bam) ]


    main:

    BOWTIE2_BUILD(metagenome) // currently waiting to see if nf-core will update to include a meta map input
    
    BOWTIE2_BUILD.out.index
      .combine(reads)
      .set{ bowtie2_align_input_ch }
    
    BOWTIE2_ALIGN(bowtie2_align_input_ch)    // currently waiting to see if nf-core will update to include a meta map input
  
    emit:
      sam = BOWTIE2_ALIGN.out.bed
}
