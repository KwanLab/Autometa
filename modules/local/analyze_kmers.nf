
process ANALYZE_KMERS {
    tag "Counting kmers for ${meta.id}"
    label 'process_medium'

    conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("*kmers.tsv")           , emit: counts
        tuple val(meta), path("*kmers.normalized.tsv"), emit: normalized
        tuple val(meta), path("*kmers.embedded.tsv")  , emit: embedded
        path 'versions.yml'                          , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        autometa-kmers \\
            --fasta ${metagenome} \\
            --kmers "${prefix}.kmers.tsv" \\
            --size "${params.kmer_size}" \\
            --norm-output "${prefix}.kmers.normalized.tsv" \\
            --norm-method "${params.norm_method}" \\
            --pca-dimensions "${params.pca_dimensions}" \\
            --embedding-output "${prefix}.kmers.embedded.tsv" \\
            --embedding-method "${params.embedding_method}" \\
            --embedding-dimensions "${params.embedding_dimensions}" \\
            --cpus "${task.cpus}" \\
            --seed 42

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
