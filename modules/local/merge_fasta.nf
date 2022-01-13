// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_FASTA {
    tag "Merging ${meta.id} FASTA"
    label 'process_low'

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

    conda (params.enable_conda ? "bioconda::seqkit=0.16.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0"
    }

    input:
        tuple val(meta), path("?")
        val extension

    output:
        tuple val(meta), path("${meta.id}.${extension}"), emit: merged
        path '*.version.txt'                            , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        # If errors occur because of issues with symlinks,
        # try:  cat * | seqkit sort -n > "${meta.id}.${extension}"
        seqkit sort -n * > "${meta.id}.${extension}"
        seqkit version | sed 's/seqkit v//g' > ${software}.version.txt
        """
}
