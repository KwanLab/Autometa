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
    kmer_coverage_results
    gc_content
    markers_results
    taxon_assignment_results
    binning_column
    

    main: 
        BINNING(kmers_embedded, kmer_coverage_results, gc_content, markers_results, taxon_assignment_results)
        UNCLUSTERED_RECRUITMENT(kmers_normalized, kmer_coverage_results, BINNING.out.binning, markers_results, taxon_assignment_results)
        BINNING_SUMMARY(BINNING.out.main, markers_results, metagenome, binning_column)
  
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