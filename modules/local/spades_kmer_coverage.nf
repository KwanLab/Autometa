process SPADES_KMER_COVERAGE {
    tag "${meta.id}"
    label 'process_low'

    conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("coverage.tsv")     , emit: coverage
        path  'versions.yml'                     , emit: versions

    when:
        meta.cov_from_assembly.equals('spades')

    script:
        """
        autometa-coverage \\
            --assembly ${metagenome} \\
            --from-spades \\
            --out "coverage.tsv"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
