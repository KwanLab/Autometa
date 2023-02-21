process BEDTOOLS_GENOMECOV {
    tag "${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "quay.io/biocontainers/bedtools:2.30.0--hc088bd4_0"
    }

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path("alignments.bed"), emit: bed
        path  "*.version.txt"                  , emit: version

    when:
        meta.cov_from_assembly.equals('0')

    when:
        task.ext.when == null || task.ext.when

    script:
        def software = getSoftwareName(task.process)
        """
        bedtools \\
            genomecov \\
            -ibam ${bam} \\
            $options.args  > alignments.bed

        bedtools --version | sed -e "s/bedtools v//g" > ${software}.version.txt
        """
}
