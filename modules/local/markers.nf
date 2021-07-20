// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

// TODO: For faster results/ les I/O this could be replaced with hmmsearch
process MARKERS {
    tag "Finding markers for ${meta.id}"
    label "process_low"

    conda (params.enable_conda ? "autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jason-c-kwan/autometa:nfcore"
    }

    input:
        tuple val(meta), path(orfs)
        //path(hmmdb) currently only inside docker
        //path(cutoffs) currently only inside docker

    output:
        tuple val(meta), path("${meta.id}.markers.tsv"), emit: markers_tsv
        path  '*.version.txt'                          , emit: version

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def software = getSoftwareName(task.process)
    if (params.enable_conda)
    """
    exit 1
    """
    else
    """
    autometa-markers \\
        --orfs $orfs \\
        --hmmscan ${meta.id}.hmmscan.tsv \\
        --out ${meta.id}.markers.tsv \\
        --kingdom ${params.kingdom} \\
        --parallel \\
        --cpus ${task.cpus} \\
        --seed 42 \\
        --hmmdb "/scratch/dbs/markers/${params.kingdom}.single_copy.hmm" \\
        --cutoffs "/scratch/dbs/markers/${params.kingdom}.single_copy.cutoffs"

    echo "TODO" > autometa.version.txt"""
}
