process PREPARE_LCA {
    tag "Preparing db cache for ${dbtype}"
    label 'process_medium'

    conda "bioconda::autometa"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        path taxdump_files // instead of passing to --dbdir, stage and pass '.'
        val dbtype

    output:
        path "cache"           , emit: cache
        path 'versions.yml'    , emit: versions

    when:
        task.ext.when == null || task.ext.when

    // storeDir = (dbtype == 'gtdb') ? params.gtdb_dir : (dbtype == 'ncbi' ? params.lca_dir : null)

    script:
        """
        # https://autometa.readthedocs.io/en/latest/scripts/taxonomy/lca.html
        autometa-taxonomy-lca \\
            --blast . \\
            --lca-output . \\
            --dbdir . \\
            --dbtype ${dbtype} \\
            --cache cache \\
            --only-prepare-cache

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
