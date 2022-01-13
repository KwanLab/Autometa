// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPADES_KMER_COVERAGE {
    tag "${meta.id}"
    label 'process_low'

    publishDir "${meta.id}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(
                filename:filename,
                options:params.options,
                publish_dir:getSoftwareName(task.process),
                meta:[:],
                publish_by_meta:[]
            )
        }

    conda (params.enable_conda ? "autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("coverage.tsv")     , emit: coverages
        path  '*.version.txt'                     , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-coverage \\
            --assembly ${metagenome} \\
            --from-spades \\
            --out "coverage.tsv"

        echo "TODO" > autometa.version.txt
        """
}
