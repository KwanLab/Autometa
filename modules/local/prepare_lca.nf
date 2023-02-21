process PREPARE_LCA {
    tag "Preparing db cache from ${blastdb_dir}"
    label 'process_medium'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    storeDir 'db/lca'
    cache 'lenient'

    input:
        path(blastdb_dir)

    output:
        path "cache"           , emit: cache
        path '*.version.txt'   , emit: version

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        autometa-taxonomy-lca \\
            --blast . \\
            --lca-output . \\
            --dbdir ${blastdb_dir} \\
            --dbtype ncbi \\
            --cache cache \\
            --only-prepare-cache
        autometa --version | sed -e "s/autometa: //g" > software.version.txt
        """
}
