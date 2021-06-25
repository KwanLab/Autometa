// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPLIT_KINGDOMS {
    tag "Splitting votes into kingdoms for ${meta.id}"
    label 'process_medium'
    
    publishDir "${params.interim_dir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/TODO"
    } else {
         container "jason-c-kwan/autometa:nfcore"
         containerOptions = "-v ${params.single_db_dir}:/ncbi:rw"
    }

    input:
      tuple val(meta), path(votes), path(assembly)

    output:
      path "${meta.id}.taxonomy.tsv", emit: taxonomy
      path "${meta.id}.bacteria.fna", emit: bacteria
      path "${meta.id}.archaea.fna", emit: archaea, optional:true

    """
    autometa-taxonomy \
      --votes ${votes} \
      --output . \
      --prefix ${meta.id} \
      --split-rank-and-write superkingdom \
      --assembly ${assembly} \
      --ncbi /ncbi
    """
}
