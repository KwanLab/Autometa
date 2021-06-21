#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.rev_reads = null
params.fwd_reads = null

params.length_table_options       = [:]
params.align_reads_options        = [:]
params.sort_reads_options         = [:]
params.genome_coverage_options    = [:]

include { LENGTH_TABLE } from './../../modules/local/length_table.nf'   addParams( options: params.length_table_options )
include { ALIGN_READS } from './../../modules/local/align_reads.nf'     addParams( options: params.align_reads_options )
include { SORT_READS} from './../../modules/local/sort_reads.nf'        addParams( options: params.sort_reads_options )
include { GENOMECOV} from './../../subworkflows/local/genome_coverage.nf'  addParams( options: params.genome_coverage_options )

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
