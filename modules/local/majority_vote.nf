
process MAJORITY_VOTE {
    tag "Performing taxon majority vote on ${meta.id}"
    label 'process_medium'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(lca)
        path taxdump_files // instead of passing to --dbdir, stage and pass '.'

    output:
        tuple val(meta), path("*votes.tsv"), emit: votes
        path  'versions.yml'             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        autometa-taxonomy-majority-vote \\
            --lca ${lca} \\
            --output ${prefix}.votes.tsv \\
            --dbdir . \\
            --dbtype ncbi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
