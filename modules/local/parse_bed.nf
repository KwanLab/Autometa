process PARSE_BED {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bedtools:2.30.0--hc088bd4_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(bed)

    output:
        tuple val(meta), path("*coverage.tsv.gz")   , emit: coverage
        path  "versions.yml"                        , emit: versions

    when:
        meta.cov_from_assembly.equals('0')
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        # NOTE: Here we supply an argument to ibam to prevent raising an error
        # However, bed is the only arg required for nextflow since bed is generated from BEDTOOLS_GENOMECOV...
        autometa-bedtools-genomecov \\
            --ibam . \\
            --bed $bed \\
            --output ${prefix}.coverage.tsv

        gzip -6  ${prefix}.coverage.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
