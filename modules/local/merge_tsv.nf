// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_TSV_WITH_HEADERS {
    tag "Merging files from parallel split for ${meta.id}"
    label 'process_low'

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::autometa" : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
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
