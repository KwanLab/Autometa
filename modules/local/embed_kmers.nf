// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process EMBED_KMERS {
    tag "PCA dims:${params.pca_dimensions}, dims:${params.embedding_dimensions}, method:${params.embedding_method}, sample:${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(normalized)

    output:
        tuple val(meta), path("kmers.embedded.tsv")  , emit: embedded
        path  '*.version.txt'                        , emit: version

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

        echo "TODO" > autometa.version.txt
        """
}
