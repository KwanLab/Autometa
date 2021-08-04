// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPADES_KMER_COVERAGE {
    tag "Calculating k-mer coverage for ${meta.id}"
    label 'process_low'

    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jason-c-kwan/autometa:${params.autometa_image}"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("${meta.id}.coverage")     , emit: coverages
        path  '*.version.txt'                            , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-coverage \\
            --assembly ${metagenome} \\
            --from-spades \\
            --out "${meta.id}.coverage"

        echo "TODO" > autometa.version.txt
        """
}
