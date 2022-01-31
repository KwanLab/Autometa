// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process NORMALIZE_KMERS {
    tag "method:${params.norm_method}, sample:${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(counts)

    output:
        tuple val(meta), path("kmers.normalized.tsv"), emit: normalized
        path  '*.version.txt'                        , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-kmers \\
            --kmers $counts \\
            --norm-output "kmers.normalized.tsv" \\
            --norm-method "${params.norm_method}" \\
            --seed 42

        echo "TODO" > autometa.version.txt
        """
}
