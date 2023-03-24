#!/usr/bin/env nextflow

process ALIGN_READS {
    tag "Aligning reads to ${meta.id}"
    label 'process_high'

    conda "bioconda::autometa=${params.autometa_image_tag}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome), path(fwd_reads), path(rev_reads)

    output:
        tuple val(meta), path("*.alignments.bam"), emit: bam
        path "*.db*.bt2"                         , emit: bt2_db
        path "versions.yml"                      , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        bowtie2-build \\
            ${args} \\
            ${metagenome} \\
            ${prefix}.db \\
            --threads ${task.cpus}

        bowtie2 \\
            -x ${meta.id}.db \\
            ${args2} \\
            -p ${task.cpus} \\
            -1 $fwd_reads \\
            -2 $rev_reads \\
                | samtools view ${args} -bS - \\
                    | samtools sort -o ${prefix}.alignments.bam


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bowtie2: \$(echo \$(bowtie2 --version 2>&1) | sed 's/^.*bowtie2-align-s version //; s/ .*\$//')
        END_VERSIONS
        """
}
