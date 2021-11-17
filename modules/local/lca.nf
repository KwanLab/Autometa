// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LCA {
    tag "Finding LCA for ${meta.id}"
    label 'process_high'
    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jason-c-kwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(blast)
        path(blastdb_dir)

    output:
        tuple val(meta), path("${meta.id}.lca.tsv"), emit: lca
        path  '*.version.txt'                      , emit: version


    script:
        def software = getSoftwareName(task.process)
        """
        autometa-taxonomy-lca --blast ${blast} --dbdir ${blastdb_dir} --output ${meta.id}.lca.tsv
        echo "TODO" > autometa.version.txt
        """
}