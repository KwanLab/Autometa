#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { AUTOMETA } from './nextflow/autometa.nf'

// Below listed parameters should be provided by the parameters.config file
// Available here: https://raw.githubusercontent.com/KwanLab/Autometa/dev/nextflow/parameters.config
// 

// Check User data inputs
params.metagenome = null
if ( !params.metagenome || params.metagenome instanceof Boolean )
error """
You must supply the `metagenome` parameter in the config or on the command line!
e.g.
nextflow run main.nf -c parameters.config --metagenome "</path/to/your/metagenome(s)>"
"""
// Where to store intermediate and final results:
params.interim = null
if ( !params.interim || params.interim instanceof Boolean )
error """
You must supply the `--interim` parameter in the config or on the command line!
e.g.
nextflow run main.nf -c parameters.config --interim "</directory/path/to/store/interimediate/results>""
"""
params.processed = null
if ( !params.processed || params.processed instanceof Boolean )
error """
You must supply the `--processed` parameter in the config or on the command line!
e.g.
nextflow run main.nf -c parameters.config --processed "</directory/path/to/final/results>""
"""


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
 -----------------------------------------------------
"""

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
