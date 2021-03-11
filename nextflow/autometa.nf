#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { LENGTH_FILTER; KMERS; COVERAGE; ORFS; MARKERS } from './common-tasks.nf'
include { TAXON_ASSIGNMENT } from './taxonomy-tasks.nf'
include { BINNING; UNCLUSTERED_RECRUITMENT } from './binning-tasks.nf'

workflow AUTOMETA {
  take:
    metagenome

  main:
    // Perform various annotations on provided metagenome
    LENGTH_FILTER(metagenome)
    COVERAGE(LENGTH_FILTER.out)
    ORFS(LENGTH_FILTER.out)
    MARKERS(ORFS.out.prots)
    // Perform taxon assignment with filtered metagenome
    TAXON_ASSIGNMENT(LENGTH_FILTER.out, ORFS.out.prots)
    // Now perform binning with all of our annotations.
    KMERS(TAXON_ASSIGNMENT.out.bacteria)
    // KMERS(TAXON_ASSIGNMENT.out.archaea) ... for case of performing binning on archaea
    BINNING(KMERS.out.normalized, COVERAGE.out, MARKERS.out, TAXON_ASSIGNMENT.out.taxonomy)
    // Then unclustered recruitment of any unclustered contigs using binning assignments from above.
    UNCLUSTERED_RECRUITMENT(KMERS.out.normalized, COVERAGE.out, BINNING.out.binning, MARKERS.out, TAXON_ASSIGNMENT.out.taxonomy)


  emit:
    binning = BINNING.out.binning
    recruitment = UNCLUSTERED_RECRUITMENT.out
    all_binning_results = BINNING.out.binning | mix(UNCLUSTERED_RECRUITMENT.out) | collect
}
