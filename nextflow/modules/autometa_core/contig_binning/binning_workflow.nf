#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { BINNING } from './process/binning.nf'
include { UNCLUSTERED_RECRUITMENT } from './process/unclustered_recruitment.nf'
include { BINNING_SUMMARY } from './process/binning_summary.nf'


workflow BIN_CONTIGS {

    take:
    metagenome
    kmers_embedded
    kmers_normalized
    coverage
    gc_content
    markers
    taxon_assignments
    binning_column
    

    main: 
        BINNING(kmers_embedded, coverage, gc_content, markers, taxon_assignments)
        UNCLUSTERED_RECRUITMENT(kmers_normalized, coverage, BINNING.out.binning, markers, taxon_assignments)
        BINNING_SUMMARY(BINNING.out.main, markers, metagenome, binning_column)
  
    emit:
        binning = BINNING.out.binning
        binning_main = BINNING.out.main
        recruitment = UNCLUSTERED_RECRUITMENT.out.binning
        recruitment_main = UNCLUSTERED_RECRUITMENT.out.main
        all_binning_results = BINNING.out.binning | mix(UNCLUSTERED_RECRUITMENT.out) | collect
        summary_stats = BINNING_SUMMARY.out.stats
        summary_taxa = BINNING_SUMMARY.out.taxonomies
        metabins = BINNING_SUMMARY.out.metabins

}
