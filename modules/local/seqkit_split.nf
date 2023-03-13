/*
=======================
TODO: Not yet implemented
=======================

process SEQKIT_SPLIT {
    tag "Splitting $meta.id for parallel processing"
    label 'process_medium'

    conda "bioconda::seqkit=0.16.1"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0"
    }

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("outfolder/*")    , emit: fasta
        path "versions.yml"                     , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        seqkit \\
            split \\
            ${fasta} \\
            ${options.args} \\
            ${options.args2} \\
            -O outfolder

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        END_VERSIONS
        """
}

*/
