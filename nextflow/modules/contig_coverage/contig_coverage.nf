#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.rev_reads = null
params.fwd_reads = null

include { LENGTH_TABLE } from './process/length_table.nf'
include { ALIGN_READS } from './process/align_reads.nf'
include { SORT_READS} from './process/sort_reads.nf'
include { GENOMECOV} from './process/genome_coverage.nf'


workflow READ_COVERAGE {
  take:
    metagenome
    fwd_reads
    rev_reads
    se_reads


  main:
    LENGTH_TABLE(metagenome)
    ALIGN_READS(metagenome, fwd_reads, rev_reads, se_reads)
    SORT_READS(ALIGN_READS.out)
    GENOMECOV(SORT_READS.out, LENGTH_TABLE.out)

  emit:
    bed = GENOMECOV.out.bed
    coverage = GENOMECOV.out.coverage
}

workflow {
  take:
    metagenome
    fwd_reads
    rev_reads
    se_reads
  
  main:
    READ_COVERAGE(metagenome, fwd_reads, rev_reads, se_reads)

}


