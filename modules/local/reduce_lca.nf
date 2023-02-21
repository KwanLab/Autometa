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
        path(blastdb_dir)
        path(lca_cache)

    output:
        tuple val(meta), path("lca.tsv"), emit: lca
        path "lca_error_taxids.tsv"     , emit: error_taxids
        path "sseqid2taxid.tsv"         , emit: sseqid_to_taxids
        path '*.version.txt'            , emit: version


    when:
        task.ext.when == null || task.ext.when

    script:
        """
        autometa-taxonomy-lca \\
            --blast ${blast} \\
            --dbdir ${blastdb_dir} \\
            --dbtype ncbi \\
            --cache ${lca_cache} \\
            --lca-error-taxids lca_error_taxids.tsv \\
            --sseqid2taxid-output sseqid2taxid.tsv \\
            --lca-output lca.tsv
        autometa --version | sed -e "s/autometa: //g" > software.version.txt
        """
}
