
process MAJORITY_VOTE {
    tag "Performing taxon majority vote on ${meta.id}"
    label 'process_medium'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(lca)
        path(ncbi_tax_dir)

    output:
        tuple val(meta), path("votes.tsv"), emit: votes
        path  '*.version.txt'             , emit: version

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        autometa-taxonomy-majority-vote \\
            --lca ${lca} \\
            --output votes.tsv \\
            --dbdir "${ncbi_tax_dir}" \\
            --dbtype ncbi

        autometa --version | sed -e "s/autometa: //g" > software.version.txt
        """
}
