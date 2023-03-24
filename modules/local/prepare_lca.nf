process PREPARE_LCA {
    tag "Preparing db cache from ${params.taxdump_tar_gz_dir}"
    label 'process_medium'

    conda "bioconda::autometa=${params.autometa_image_tag}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    cache 'lenient'

    input:
        path taxdump_files // instead of passing to --dbdir, stage and pass '.'

    output:
        tuple path("./cache/level.pkl.gz"), path("./cache/occurrence.pkl.gz"), path("./cache/precomputed_lcas.pkl.gz"), path("./cache/tour.pkl.gz"), emit: pickles
        path 'versions.yml' , emit: versions

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
