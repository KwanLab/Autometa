#!/usr/bin/env nextflow


process ALIGN_READS {
    tag "Aligning reads to ${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

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

    script:
        """
        bowtie2-build \\
            ${options.args} \\
            ${metagenome} \\
            ${meta.id}.db

        bowtie2 \\
            -x ${meta.id}.db \\
            ${options.args2} \\
            -p ${task.cpus} \\
            -S alignments.sam \\
            -1 $fwd_reads \\
            -2 $rev_reads

        echo \$(bowtie2 --version 2>&1) | sed -n 's/^.*bowtie2-align-s version //p; s/ .*\$//' > bowtie2.version.txt
        """
}
