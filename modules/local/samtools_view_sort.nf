// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SAMTOOLS_VIEW_AND_SORT {
    tag "$meta.id"
    label 'process_medium'

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::samtools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0"
    } else {
        container "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("*.bam"), emit: bam
        path "*.version.txt"          , emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        """
        samtools view ${options.args} -@ ${task.cpus} ${bam} |\
            samtools sort ${options.args2} -@ ${task.cpus} -o ${prefix}.bam -T $prefix

        echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
        """
}
