#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { DIAMOND } from './process/diamond.nf'

workflow DIAMOND_BLASTP {
    take:
      orfs

    main:
      DIAMOND(orfs)
    emit:
      blastp_table = DIAMOND.out.blastp_table // output '${orfs.simpleName}.blastp.tsv'; BLAST fmt 6 table
}

