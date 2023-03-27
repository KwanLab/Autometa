// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process MAJORITY_VOTE {
    tag "Performing taxon majority vote on ${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(lca)
        path(ncbi_tax_dir)

    output:
        tuple val(meta), path("votes.tsv"), emit: votes
        path  '*.version.txt'             , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-taxonomy-majority-vote \\
            --lca ${lca} \\
            --output votes.tsv \\
            --dbdir "${ncbi_tax_dir}" \\
            --dbtype ncbi

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
