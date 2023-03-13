/*
=======================
TODO: Not yet implemented
=======================

process MERGE_FASTA {
    tag "Merging ${meta.id} FASTA"
    label 'process_low'

    conda "bioconda::seqkit=0.16.1"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0"
    }

    input:
        tuple val(meta), path("?")
        val extension

    output:
        tuple val(meta), path("${meta.id}.${extension}"), emit: merged
        path 'versions.yml'                             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        # If errors occur because of issues with symlinks,
        # try:  cat * | seqkit sort -n > "${meta.id}.${extension}"
        seqkit sort -n * > "${meta.id}.${extension}"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        END_VERSIONS
        """
}

*/
