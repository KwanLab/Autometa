#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Input
params.kingdom = "bacteria"
params.length_cutoff = 3000
params.kmer_size = 4
params.completeness = 20.0
params.purity = 90.0
params.metagenome = null
if ( !params.metagenome || params.metagenome instanceof Boolean ) error "</path(s)/to/metagenome.fna>"

// Where to store intermediate and final results:
params.interim = null
if ( !params.interim || params.interim instanceof Boolean ) error "/path/to/directory/to/store/user/interimediate/results>"

params.processed = null
if ( !params.processed || params.processed instanceof Boolean ) error "</path/to/directory/to/store/user/final/results>"

// Databases
params.ncbi_database = "$HOME/Autometa/autometa/databases/ncbi"
params.diamond_database = "$HOME/Autometa/autometa/databases/ncbi/nr.dmnd"
params.markers_database = "$HOME/Autometa/autometa/databases/markers"
// Additional runtime settings
params.cpus = 2


log.info """

 Autometa - Automated Extraction of Genomes from Shotgun Metagenomes
 =====================================================
 projectDir                         : ${workflow.projectDir}
 -----------------------------------------------------
 Data
 -----------------------------------------------------
 metagenome                         : ${params.metagenome}
 interim                            : ${params.interim}
 processed                          : ${params.processed}
 -----------------------------------------------------
 Parameters
 -----------------------------------------------------
 cpus                               : ${params.cpus}
 length_cutoff                      : ${params.length_cutoff}
 kmer_size                          : ${params.kmer_size}
 kmer_norm_method                   : ${params.kmer_norm_method}
 kmer_pca_dimensions                : ${params.kmer_pca_dimensions}
 kmer_embed_method                  : ${params.kmer_embed_method}
 kmer_embed_dimensions              : ${params.kmer_embed_dimensions}
 clustering_method                  : ${params.clustering_method}
 classification_kmer_pca_dimensions : ${params.classification_kmer_pca_dimensions}
 classification_method              : ${params.classification_method}
 completeness                       : ${params.completeness}
 purity                             : ${params.purity}
 kingdom                            : ${params.kingdom}
 -----------------------------------------------------
 Databases
 -----------------------------------------------------
 ncbi_database                      : ${params.ncbi_database}
 diamond_database                   : ${params.diamond_database}
 markers_database                   : ${params.markers_database}
 -----------------------------------------------------
"""

// Note: It is required to include these processes after parameter definitions
// so they take on the provided values
include { LENGTH_FILTER; KMERS; COVERAGE; ORFS; MARKERS } from './nextflow/common-tasks.nf'
include { TAXON_ASSIGNMENT } from './nextflow/taxonomy-tasks.nf'
include { BINNING; UNCLUSTERED_RECRUITMENT } from './nextflow/binning-tasks.nf'

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

workflow {
  Channel
    .fromPath(params.metagenome, checkIfExists: true, type: 'file')
    .set{unfiltered_metagenome_ch}

  AUTOMETA(unfiltered_metagenome_ch)
}

/*
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
