// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PREPARE_LCA {
    tag "Preparing db cache for autometa-taxonomy-lca"
    label 'process_high'
    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jason-c-kwan/autometa:${params.autometa_image_tag}"
    }

    cache 'lenient'

    input:
        path(blastdb_dir)

    output:
        path "lca_cache"       , emit: cache
        path '*.version.txt'   , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-taxonomy-lca \\
            --blast . \\
            --lca-output . \\
            --dbdir ${blastdb_dir} \\
            --cache lca_cache \\
            --only-prepare-cache
        echo "TODO" > autometa.version.txt
        """
}