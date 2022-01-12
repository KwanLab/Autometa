// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPLIT_KINGDOMS {
    tag "Splitting votes into kingdoms for ${meta.id}"
    label 'process_medium'

    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(assembly), path(votes)
        path(ncbi_tax_dir)

    output:
        tuple val(meta), path("${meta.id}.taxonomy.tsv"), emit: taxonomy
        tuple val(meta), path("${meta.id}.bacteria.fna"), emit: bacteria, optional: true
        tuple val(meta), path("${meta.id}.archaea.fna") , emit: archaea, optional: true
        path  '*.version.txt'                           , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-taxonomy \\
            --votes "${votes}" \\
            --output . \\
            --prefix "${meta.id}" \\
            --split-rank-and-write superkingdom \\
            --assembly "${assembly}" \\
            --ncbi "${ncbi_tax_dir}"

        echo "TODO" > autometa.version.txt
        """
}
