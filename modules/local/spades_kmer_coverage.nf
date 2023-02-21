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
        path  '*.version.txt'                     , emit: version

    when:
        meta.cov_from_assembly.equals('spades')
        task.ext.when == null || task.ext.when

    script:
        """
        autometa-coverage \\
            --assembly ${metagenome} \\
            --from-spades \\
            --out "coverage.tsv"

        autometa --version | sed -e "s/autometa: //g" > software.version.txt
        """
}
