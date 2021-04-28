#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.interim_dir = "</path/to/store/user/interimediate/results>"
params.outdir = "</path/to/store/user/final/results>"
params.ncbi_database = "$HOME/Autometa/autometa/databases/ncbi"
params.cpus = 2

include { DIAMOND } from './process/diamond.nf'
include { LCA } from './process/lca.nf'
include { MAJORITY_VOTE } from './process/majority_vote.nf'
include { SPLIT_KINGDOMS } from './process/split_kingdoms.nf'


// Autometa taxon assignment workflow
workflow TAXON_ASSIGNMENT {
    take:
      assembly
      orfs

    main:
      DIAMOND(orfs)
      LCA(DIAMOND.out)
      MAJORITY_VOTE(LCA.out)
      SPLIT_KINGDOMS(MAJORITY_VOTE.out, assembly)

    emit:
      taxonomy = SPLIT_KINGDOMS.out.taxonomy
      bacteria = SPLIT_KINGDOMS.out.bacteria
      archaea = SPLIT_KINGDOMS.out.archaea
      orf_votes = LCA.out
      contig_votes = MAJORITY_VOTE.out
}
