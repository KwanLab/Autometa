process PARSE_BED {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(bed)

    output:
        tuple val(meta), path("coverage.tsv"), emit: coverage
        path  "*.version.txt"                , emit: version

    when:
        meta.cov_from_assembly.equals('0')
        task.ext.when == null || task.ext.when

    script:
        def software = getSoftwareName(task.process)
        """
        # NOTE: Here we supply an argument to ibam to prevent raising an error
        # However, bed is the only arg required for nextflow since bed is generated from BEDTOOLS_GENOMECOV...
        autometa-bedtools-genomecov \\
            --ibam . \\
            --bed $bed \\
            --output coverage.tsv

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
