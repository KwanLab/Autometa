process EMBED_KMERS {
    tag "PCA dims:${params.pca_dimensions}, dims:${params.embedding_dimensions}, method:${params.embedding_method}, sample:${meta.id}, taxon:${meta.taxon}"
    label 'process_medium'

    conda "bioconda::autometa=${params.autometa_image_tag}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }
    input:
        tuple val(meta), path(normalized)

    output:
        // optional because autometa will exit with 0 if not enough contigs are present
        // don't want to raise an error because even if ignored users will see a Note and be confused
        tuple val(meta), path("*kmers.embedded.tsv.gz") , emit: embedded, optional:true
        path  'versions.yml'                            , emit: versions

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

        # gzip if file exists
        find ./ -name "${prefix}.kmers.embedded.tsv" -type f -exec gzip -6 {} +

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}

