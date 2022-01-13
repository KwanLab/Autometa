// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ANALYZE_KMERS {
    tag "Counting kmers for ${meta.id}"
    label 'process_medium'
    publishDir "${meta.id}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(
                filename:filename,
                options:params.options,
                publish_dir:getSoftwareName(task.process),
                meta:[:],
                publish_by_meta:[]
            )
        }

    conda (params.enable_conda ? "autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("kmers.tsv")           , emit: counts
        tuple val(meta), path("kmers.normalized.tsv"), emit: normalized
        tuple val(meta), path("kmers.embedded.tsv")  , emit: embedded
        path  '*.version.txt'                        , emit: version

    script:
        // Add soft-links to original FastQs for consistent naming in pipeline
        def software = getSoftwareName(task.process)
        """
        autometa-kmers \\
            --fasta ${metagenome} \\
            --kmers "kmers.tsv" \\
            --size "${params.kmer_size}" \\
            --norm-output "kmers.normalized.tsv" \\
            --norm-method "${params.norm_method}" \\
            --pca-dimensions "${params.pca_dimensions}" \\
            --embedding-output "kmers.embedded.tsv" \\
            --embedding-method "${params.embedding_method}" \\
            --embedding-dimensions "${params.embedding_dimensions}" \\
            --cpus "${task.cpus}" \\
            --seed 42

        echo "TODO" > autometa.version.txt
        """
}
