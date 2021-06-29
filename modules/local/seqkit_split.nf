// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)
params.num_splits = 20
process SEQKIT_SPLIT {
    tag "$meta.id"
    label 'process_medium'
    
    // no publishdir

    conda (params.enable_conda ? "bioconda::seqkit=0.16.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0"
    }

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("outfolder/*")    , emit: fasta
    path "*.version.txt"                    , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
 
    """
    seqkit \\
        split \\
        ${fasta} \\
        -p ${params.num_splits} \\
        ${options.args} \\
        -O outfolder 

    seqkit version | sed 's/seqkit v//g' > ${software}.version.txt
    """
}
