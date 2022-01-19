// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PARSE_BED {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(bam), path(lengths), path(bed_out)

    output:
        tuple val(meta), path("${meta.id}.coverage.tsv"), emit: coverage
        path  "*.version.txt"                           , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-bedtools-genomecov \\
            --ibam $bam \\
            --lengths $lengths \\
            --bed $bed_out \\
            --output ${meta.id}.coverage.tsv

        echo "TODO" > autometa.version.txt
        """
}
