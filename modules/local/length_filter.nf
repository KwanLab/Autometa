// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LENGTH_FILTER {
    tag "Removing ${meta.id} contigs less than ${params.length_cutoff}"
    label 'process_medium'
   
    publishDir "${params.interim_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
    
    conda (params.enable_conda ? "autometa" : null)
     if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
     } else {
         container "jason-c-kwan/autometa:nfcore"
     }

    input:
        tuple val(meta), path(metagenome)
  
    output:
      tuple val(meta), path("${meta.id}.filtered.fna"), emit: fasta
      tuple val(meta), path("${meta.id}.stats.tsv"), emit: stats
      tuple val(meta), path("${meta.id}.gc_content.tsv"), emit: gc_content
  
    script:
    def software = getSoftwareName(task.process)
    // python code fails if input file ends with .filtered.fna 
    if ( "${meta.id}" == "${meta.id}.filtered.fna" )
    """
    autometa-length-filter \
      --assembly $metagenome \
      --cutoff ${params.length_cutoff} \
      --output-fasta "${meta.id}.autometa_filtered.fna" \
      --output-stats ${meta.id}.stats.tsv \
      --output-gc-content ${meta.id}.gc_content.tsv
    """
    else 
    """
    autometa-length-filter \
      --assembly $metagenome \
      --cutoff ${params.length_cutoff} \
      --output-fasta "${meta.id}.filtered.fna" \
      --output-stats ${meta.id}.stats.tsv \
      --output-gc-content ${meta.id}.gc_content.tsv
  
    """
}
