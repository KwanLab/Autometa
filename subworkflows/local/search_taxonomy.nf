#!/usr/bin/env nextflow

nextflow.enable.dsl=2


// This workflow is used instead DIAMOND_BLASTP directly to provide flexibility
// in possibly using other programs in the future 

include { DIAMOND as DIAMOND_BLASTP } from './../../modules/local/fix_diamond.nf'

workflow SEARCH_TAXONOMY {
    take:
      orfs

    main:
      DIAMOND_BLASTP(orfs)
  
    emit:
      blastp_table = DIAMOND_BLASTP.out.blastp_table // output '${orfs.simpleName}.blastp.tsv'; BLAST fmt 6 table
}
