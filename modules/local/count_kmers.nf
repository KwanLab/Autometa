
process COUNT_KMERS {
    tag "Counting ${params.kmer_size}-mers for ${meta.id}"
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
        tuple val(meta), path("kmers.tsv")           , emit: counts
        path  'versions.yml'                        , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        autometa-kmers \\
            --fasta $metagenome \\
            --kmers "kmers.tsv" \\
            --size "${params.kmer_size}" \\
            --cpus "${task.cpus}" \\
            --seed 42

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
