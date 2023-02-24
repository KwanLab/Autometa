process REDUCE_LCA {
    tag "Finding LCA for ${meta.id}"
    label 'process_medium'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(blast)
        path taxdump_files // instead of passing to --dbdir, stage and pass '.'
        path lca_cache
        path prot_accession2taxid

    output:
        tuple val(meta), path("lca.tsv"), emit: lca
        path "lca_error_taxids.tsv"     , emit: error_taxids
        path "sseqid2taxid.tsv"         , emit: sseqid_to_taxids
        path 'versions.yml'            , emit: versions


    when:
        task.ext.when == null || task.ext.when

    script:
        """
        autometa-taxonomy-lca \\
            --blast ${blast} \\
            --dbdir . \\
            --dbtype ncbi \\
            --cache ${lca_cache} \\
            --lca-error-taxids lca_error_taxids.tsv \\
            --sseqid2taxid-output sseqid2taxid.tsv \\
            --lca-output lca.tsv
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
