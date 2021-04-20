#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
                         Autometa   
========================================================================================
 Autometa's Nextflow Analysis Pipeline
 #### Homepage
 https://github.com/KwanLab/Autometa
 #### Documentation
 https://autometa.readthedocs.io/en/latest/
----------------------------------------------------------------------------------------
*/

// Below listed parameters should be provided by the parameters.config file
// Available here: https://raw.githubusercontent.com/KwanLab/Autometa/dev/nextflow/parameters.config
//

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////
def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run Autometa --input '*_R{1,2}.fastq.gz' -profile docker

    Note: Pay attention to the use of "--" vs "-" in the below arguments.

    Mandatory arguments:
      --input_fna  [file]            path to fasta.gz containing assembled contigs
      --interim_dir [file]           path to directory where intermediate (temporary) 
                                       files will be written/read
      --outdir [file]                The output directory where the results will be saved
      --result_mode [str]           Mode for publishing results in the output directory. 
                                       Available: symlink, rellink, link, copy, copyNoFollow,
                                       move (Default: copy)

    Other Autometa arguments:
    Metagenome Length filtering
      --length_cutoff [int]             Minimum contig length to use as input to Autometa
    Kmer counting/normalization/embedding
      --kmer_size [int]                 default: 5
      --kmer_norm_method [str]          default: "am_clr"
                                        choices: "am_clr", "clr", "ilr"
      --kmer_pca_dimensions [int]       default: 2
      --kmer_embed_method [str]         default: "bhsne"
                                        choices: "sksne", "bhsne", "umap"
      --kmer_embed_dimensions [int]     default: 2
    Binning parameters
      --kingdom  [str]                  default: "bacteria"
      --binning_starting_rank           default: "superkingdom"
                                        choices: "superkingdom", "phylum", "class", "order", "family", "genus", "species"
      --clustering_method [str]         default: "dbscan"
                                        choices: "dbscan", "hdbscan"
      --classification_kmer_pca_dimensions [int] e.g. 50
      --classification_method           default: "decision_tree"
                                        options: "decision_tree", "random_forest"
      --completeness [float]            default: 20.0
      --purity [float]                  default: 90.0

    Optional Nextflow arguments:
      -profile [str]                  Configuration profile to use.
                                        Available: "conda", "docker", "test"
      --max_memory [str]              e.g. "128.GB"
      --max_cpus [int]                e.g. 10

    Other options:
      --email [email]                 Set this parameter to your e-mail address to get a summary
                                        e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will 
                                        automatically generate a random mnemonic

    AWSBatch options (placeholder-not yet implemented/tested):
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
      
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Path of fna files used as input:
params.input_fna = null
if ( !params.input_fna || params.input_fna instanceof Boolean )
error """
You must supply the `input_fna` parameter in the config or on the command line!
For a single input use the path
    nextflow run main.nf -c parameters.config --input_fna "</path/to/your/metagenome.fna>"
For multiple input files, glob patterns may be used:
    nextflow run main.nf -c parameters.config --input_fna "</path/to/your/metagenome_*.fna>"
"""

// Where to store final results:
params.outdir = null
if ( !params.outdir || params.outdir instanceof Boolean )
error """
You must supply the `--outdir` parameter in the config or on the command line!
e.g.
nextflow run main.nf -c parameters.config --outdir "</directory/path/to/final/results>"
"""

// Where to store intermediate results:
params.interim_dir = null
if ( !params.interim_dir )
error """
You must supply the `--interim_dir` parameter in the config or on the command line!
e.g.
nextflow run main.nf -c parameters.config --interim_dir "</directory/path/to/store/interimediate/results>"
"""
if ( params.interim_dir = true ) {
    params.interim_dir = "${params.outdir}/interim_outputs"
    println "Intermediate results will be saved to:\n${params.outdir}/interim_outputs"
}
if ( params.interim_dir = false ) {
    // TODO: Don't publish/create the copies of intermediate results
}


log.info """

 Autometa - Automated Extraction of Genomes from Shotgun Metagenomes
 =====================================================
 projectDir                         : ${workflow.projectDir}
 -----------------------------------------------------
 Data
 -----------------------------------------------------
 metagenome                         : ${params.input_fna}
 interim                            : ${params.interim_dir}
 processed                          : ${params.outdir}
 -----------------------------------------------------
 Parameters
 -----------------------------------------------------
 cpus                               : ${params.cpus}
 length_cutoff                      : ${params.length_cutoff}
 kmer_size                          : ${params.kmer_size}
 norm_method                        : ${params.norm_method}
 pca_dimensions                     : ${params.pca_dimensions}
 embedding_method                   : ${params.embedding_method}
 embedding_dimensions               : ${params.embedding_dimensions}
 clustering_method                  : ${params.clustering_method}
 classification_kmer_pca_dimensions : ${params.classification_kmer_pca_dimensions}
 classification_method              : ${params.classification_method}
 completeness                       : ${params.completeness}
 purity                             : ${params.purity}
 gc_stddev_limit                    : ${params.gc_stddev_limit}
 cov_stddev_limit                   : ${params.cov_stddev_limit}
 kingdom                            : ${params.kingdom}
 -----------------------------------------------------
 Databases
 -----------------------------------------------------
 ncbi_database                      : ${params.ncbi_database}
 -----------------------------------------------------
"""

/*
 * Retrieve the main 'AUTOMETA' worflow
 */
include { AUTOMETA } from './nextflow/autometa.nf'


workflow {
  Channel
    .fromPath(params.input_fna, checkIfExists: true, type: 'file')
    .set{unfiltered_metagenome_ch}

  AUTOMETA(unfiltered_metagenome_ch)
}

/*
 * completion handler
 */
workflow.onComplete {
	log.info ( workflow.success ? "\nDone!\n" : "Oops .. something went wrong" )
}
