process NORMALIZE_KMERS {
    tag "method:${params.norm_method}, sample:${meta.id}"
    label 'process_medium'

    conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(counts)

    output:
        tuple val(meta), path("*kmers.normalized.tsv"), emit: normalized
        path  'versions.yml'                          , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        autometa-kmers \\
            --kmers $counts \\
            --norm-output ${prefix}.kmers.normalized.tsv \\
            --norm-method ${params.norm_method} \\
            --seed 42

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
