process EMBED_KMERS {
    tag "PCA dims:${params.pca_dimensions}, dims:${params.embedding_dimensions}, method:${params.embedding_method}, sample:${meta.id}"
    label 'process_medium'

    conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    // Not enough contigs to perform embedding with current parameter settings...
    errorStrategy { task.exitStatus in 153 ? 'ignore' : 'terminate' }

    input:
        tuple val(meta), path(normalized)

    output:
        tuple val(meta), path("kmers.embedded.tsv")  , emit: embedded
        path  '*.version.txt'                        , emit: version

    when:
        task.ext.when == null || task.ext.when

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-kmers \\
            --norm-output $normalized \\
            --pca-dimensions "${params.pca_dimensions}" \\
            --embedding-output "kmers.embedded.tsv" \\
            --embedding-method "${params.embedding_method}" \\
            --embedding-dimensions "${params.embedding_dimensions}" \\
            --cpus "${task.cpus}" \\
            --seed 42

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
