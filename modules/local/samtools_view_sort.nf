// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

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
        tuple val(meta), path(sam)

    output:
        tuple val(meta), path("alignments.bam"), emit: bam
        path "*.version.txt"                   , emit: version

    when:
        meta.cov_from_assembly.equals('0')

    script:
        def software = getSoftwareName(task.process)
        """
        samtools view ${options.args} -@ ${task.cpus} -bS ${sam} \\
            | samtools sort ${options.args2} -@ ${task.cpus} -o alignments.bam

        echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
        """
}
