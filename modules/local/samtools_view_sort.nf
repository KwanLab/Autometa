process SAMTOOLS_VIEW_AND_SORT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.13"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0"
    } else {
        container "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }

    input:
        tuple val(meta), path(sam)

    output:
        tuple val(meta), path("*.alignments.bam"), emit: bam
        path "versions.yml"                   , emit: versions

    when:
        meta.cov_from_assembly.equals('0')
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"

        """
        samtools view ${args} -@ ${task.cpus} -bS ${sam} \\
            | samtools sort ${args2} -@ ${task.cpus} -o ${prefix}.alignments.bam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
}
