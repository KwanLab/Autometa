process PREPARE_LCA {
    //tag "Preparing db cache from ${blastdb_dir}"
    label 'process_medium'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    cache 'lenient'

    input:
        path taxdump_files // instead of passing to --dbdir, stage and pass '.'

    output:
        path "cache"           , emit: cache
        path 'versions.yml'   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        # https://autometa.readthedocs.io/en/latest/scripts/taxonomy/lca.html
        autometa-taxonomy-lca \\
            --blast . \\
            --lca-output . \\
            --dbdir . \\
            --dbtype ncbi \\
            --cache cache \\
            --only-prepare-cache
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
