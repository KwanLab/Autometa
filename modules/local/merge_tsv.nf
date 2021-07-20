// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_TSV_WITH_HEADERS {
    tag "Merging files from parallel split for ${meta.id}"
    label 'process_low'
    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::autometa" : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
         container "jason-c-kwan/autometa:${params.autometa_image}"
    }

    input:
        tuple val(meta), path("?.tsv")
        val extension

    output:
    tuple val(meta), path("${meta.id}.${extension}"), emit: merged_tsv


    script:
    def software = getSoftwareName(task.process)


    """
      awk 'FNR==1 && NR!=1{next;}{print}' *.tsv > "${meta.id}.${extension}"
    """
}
