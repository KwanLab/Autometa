#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MAJORITY_VOTE {
  label 'process_medium'
  
  tag "Performing taxon majority vote on ${lca.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:rw"
  publishDir params.interim_dir, pattern: "${lca.simpleName}.votes.tsv"

  input:
    path lca

  output:
    path "${lca.simpleName}.votes.tsv"

  """
  autometa-taxonomy-majority-vote --lca ${lca} --output ${lca.simpleName}.votes.tsv --dbdir /ncbi
  """
}


// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LCA {
    tag "Finding LCA for ${meta.id}"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::autometa" : null)
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
         container "jason-c-kwan/autometa:nfcore"
    }

    input:
        tuple val(meta), path(blast)

    output:
    tuple val(meta), path("${meta.id}.lca.tsv"), emit: lca


    script:
    def software = getSoftwareName(task.process)
    
    """
    autometa-taxonomy-majority-vote \\
        --lca ${blast} \\
        --output ${meta.id}.votes.tsv \\
        --dbdir /ncbi
    """
}
