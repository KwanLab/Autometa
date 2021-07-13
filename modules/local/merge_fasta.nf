// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_FASTA {
    tag "Merging ${meta.id} FASTA"
    label 'process_low'
    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
        tuple val(meta), path("?")
        val extension

    output:
    tuple val(meta), path("${meta.id}.${extension}"), emit: merged


    script:
    def software = getSoftwareName(task.process)
    

    """
     cat * > "${meta.id}.${extension}"
    """
}
