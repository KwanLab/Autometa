#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LCA } from './process/lca.nf'
include { MAJORITY_VOTE } from './process/majority_vote.nf'
include { SPLIT_KINGDOMS } from './process/split_kingdoms.nf'


// Autometa taxon assignment workflow
workflow TAXON_ASSIGNMENT {
    take:
      metagenome
      blastp_table

    main:
      LCA(blastp_table) // output '${blast.simpleName}.lca.tsv'
      MAJORITY_VOTE(LCA.out) //output ${lca.simpleName}.votes.tsv
      SPLIT_KINGDOMS(MAJORITY_VOTE.out, metagenome) // output "${assembly.simpleName}.taxonomy.tsv" "${assembly.simpleName}.bacteria.fna" , "${assembly.simpleName}.archaea.fna"

    emit:
      taxonomy = SPLIT_KINGDOMS.out.taxonomy
      bacteria = SPLIT_KINGDOMS.out.bacteria
      archaea = SPLIT_KINGDOMS.out.archaea
      orf_votes = LCA.out
      contig_votes = MAJORITY_VOTE.out
}
