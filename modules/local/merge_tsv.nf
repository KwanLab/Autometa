/*
=======================
TODO: Not yet implemented
=======================

process MERGE_TSV_WITH_HEADERS {
    tag "Merging files from parallel split for ${meta.id}"
    label 'process_low'

    conda "bioconda::autometa"

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


    when:
        task.ext.when == null || task.ext.when

    script:
        """
        awk 'FNR==1 && NR!=1{next;}{print}' *.tsv > "${meta.id}.${extension}"
        """
}


*/

