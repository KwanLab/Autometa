// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

// nf-core version doesn't have the easy ability to map the -g flag

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
        container "jason-c-kwan/autometa:${params.autometa_image}"
    }

    input:
    tuple val(meta), path(bam), path(lengths), path(bed_out)

    output:
    tuple val(meta), path("${meta.id}.coverage.tsv"), emit: coverage
    path  "*.version.txt"                           , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    autometa-parse-bed \
        --ibam $bam \
        --lengths $lengths \
        --bed $bed_out \
        --output ${meta.id}.coverage.tsv

    echo "TODO" > autometa.version.txt
    """
}
