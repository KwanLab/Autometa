#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process GENOMECOV {
  tag "Computing genome coverage for ${bam.simpleName}"

  input:
    path bam
    path lengths

  output:
    path "${bam.simpleName}.bed.tsv", emit: bed
    path "${bam.simpleName}.coverage.tsv", emit: coverage

  """
  bedtools genomecov -ibam $bam -g $lengths > ${bam.simpleName}.bed.tsv
  autometa-parse-bed \
    --ibam $bam \
    --lengths $lengths \
    --bed ${bam.simpleName}.bed.tsv \
    --output ${bam.simpleName}.coverage.tsv
  """
}

