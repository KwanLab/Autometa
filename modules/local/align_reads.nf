#!/usr/bin/env nextflow


process ALIGN_READS {
    tag "Aligning reads to ${meta.id}"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome), path(fwd_reads), path(rev_reads)

    output:
        tuple val(meta), path("alignments.sam"), emit: sam
        path "*.db*.bt2"                       , emit: bt2_db
        path "*.version.txt"                   , emit: version

    when:
        meta.cov_from_assembly.equals('0')
        task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
        """
        bowtie2-build \\
            ${args} \\
            ${metagenome} \\
            ${meta.id}.db

        bowtie2 \\
            -x ${meta.id}.db \\
            ${args2} \\
            -p ${task.cpus} \\
            -S alignments.sam \\
            -1 $fwd_reads \\
            -2 $rev_reads

        echo \$(bowtie2 --version 2>&1) | sed -n 's/^.*bowtie2-align-s version //p; s/ .*\$//' > bowtie2.version.txt
        """
}
