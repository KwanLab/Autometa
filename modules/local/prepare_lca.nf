// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PREPARE_LCA {
    tag "Preparing db cache from ${blastdb_dir}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
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

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-taxonomy-lca \\
            --blast . \\
            --lca-output . \\
            --dbdir ${blastdb_dir} \\
            --cache cache \\
            --only-prepare-cache
        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
