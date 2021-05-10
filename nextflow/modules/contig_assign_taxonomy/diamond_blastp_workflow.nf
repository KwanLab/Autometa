#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.interim_dir = "</path/to/store/user/interimediate/results>"
params.outdir = "</path/to/store/user/final/results>"
params.ncbi_database = "$HOME/Autometa/autometa/databases/ncbi"
params.cpus = 2

include { DIAMOND } from './process/diamond.nf'

workflow DIAMOND_BLASTP {
    take:
      orfs

    main:
      DIAMOND(orfs)
    emit:
      diamond_blastp_table = DIAMOND.diamond_blastp_table // output '${orfs.simpleName}.blastp.tsv'; BLAST fmt 6 table
}
