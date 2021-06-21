// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)



process BINNING_SUMMARY {
    tag "Gathering binning summary for ${meta.id}"
    label 'low'
   
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
   
    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jason-c-kwan/autometa:nfcore"
    }
    
    input:
        tuple val(meta), path(binning_main), path(markers), path(metagenome)
        val(binning_column)

    output:
        tuple val(meta), path("${meta.id}_metabin_stats.tsv"), emit: stats
        tuple val(meta), path("${meta.id}_metabins"), emit: metabins
        tuple val(meta), path("${meta.id}_metabin_taxonomy.tsv"), emit: taxonomies, optional: true

    script:
    """
    mkdir -p ${meta.id}

    autometa-binning-summary \
      --ncbi /ncbi \
      --binning-main $binning_main \
      --markers $markers \
      --metagenome $metagenome \
      --binning-column $binning_column \
      --output-stats "${meta.id}_metabin_stats.tsv" \
      --output-taxonomy "${meta.id}_metabin_taxonomy.tsv" \
      --output-metabins "${meta.id}_metabins"
    """
}
