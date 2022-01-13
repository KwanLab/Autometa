
/*
=======================
TODO: Not yet implemented
"Cutoffs" need to be downloaded/provided to this process
=======================
*/


// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HMMER_HMMSEARCH_FILTER {
    tag "Filtering marker hmms in ${meta.id}"
    label 'process_medium'

    // if ( params.num_splits < 2 ) {
    // if running in parallel, the results are published from the process
    // that merges the individual results from this process
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
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(domtblout), path(fasta)
        path("bacteria.single_copy.cutoffs") // FIXME: This should take cutoffs corresponding to domtblout

    output:
        tuple val(meta), path("markers.tsv"), emit: markers_tsv
        path "*.version.txt"                           , emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        """
        autometa-hmmsearch-filter \\
            --domtblout "$domtblout" \\
            --cutoffs  TODO:"Cutoffs" need to be downloaded/provided to this process \\
            --seqdb "$fasta" \\
            --out "markers.tsv"

        echo "TODO" > autometa.version.txt
        """
}
