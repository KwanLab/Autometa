// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

//TODO REMOVE FILE, this has been replaced by hmmsearch
process MARKERS {
    tag "Finding markers for ${meta.id}"
    label "process_low"
    
    // copying orfs via stageInMode is required to run hmmscan (does not handle symlinks)
    //stageInMode 'copy'
   
    conda (params.enable_conda ? "autometa" : null)
     if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
     } else {
         container "jason-c-kwan/autometa:nfcore"
     }  
   
    input:
        tuple val(meta), path(orfs)
  
    output:
        tuple val(meta), path("${meta.id}.markers.tsv"), emit: markers_tsv
  
    """
    autometa-markers \\
        --orfs $orfs \\
        --hmmscan ${meta.id}.hmmscan.tsv \\
        --out ${meta.id}.markers.tsv \\
        --kingdom ${params.kingdom} \\
        --parallel \\
        --cpus ${task.cpus} \\
        --seed 42
    """
}
