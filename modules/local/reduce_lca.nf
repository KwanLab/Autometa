// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process REDUCE_LCA {
    tag "Finding LCA for ${meta.id}"
    label 'process_high'
    publishDir "${meta.id}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::autometa" : null)
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


    script:
        def software = getSoftwareName(task.process)
        """
        autometa-taxonomy-lca \\
            --blast ${blast} \\
            --dbdir ${blastdb_dir} \\
            --cache ${lca_cache} \\
            --lca-error-taxids lca_error_taxids.tsv \\
            --sseqid2taxid-output sseqid2taxid.tsv \\
            --lca-output lca.tsv
        echo "TODO" > autometa.version.txt
        """
}
