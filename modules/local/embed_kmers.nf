process EMBED_KMERS {
    tag "PCA dims:${params.pca_dimensions}, dims:${params.embedding_dimensions}, method:${params.embedding_method}, sample:${meta.id}, taxon:${meta.taxon}"
    label 'process_medium'
    // Not enough contigs to perform embedding with current parameter settings...
    errorStrategy { task.exitStatus in 8 ? 'ignore' : 'terminate' }

    conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }
    input:
        tuple val(meta), path(normalized)

    output:
        tuple val(meta), path("*kmers.embedded.tsv.gz")  , emit: embedded
        path  'versions.yml'                          , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}.${meta.taxon}"
        """
        autometa-kmers \\
            --norm-output ${normalized} \\
            --pca-dimensions ${params.pca_dimensions} \\
            --embedding-output ${prefix}.kmers.embedded.tsv \\
            --embedding-method ${params.embedding_method} \\
            --embedding-dimensions ${params.embedding_dimensions} \\
            --cpus ${task.cpus} \\
            --seed 42

        gzip -6  "${prefix}.kmers.embedded.tsv"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}

