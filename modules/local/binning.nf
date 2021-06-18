// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BINNING {
    tag "Performing Autometa binning on ${meta.id}"
    label 'low'
   
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
   
    conda (params.enable_conda ? "autometa" : null)
     if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
     } else {
         container "jason-c-kwan/autometa:nfcore"
     }
    
    input:
        tuple val(meta), path(kmers), path(coverage), path(gc_content), path(markers)
        val(taxonomy)
  
    output:
        tuple val(meta), path("${coverage.simpleName}.${params.kingdom}.binning.tsv.gz"), emit: binning
        tuple val(meta), path("${coverage.simpleName}.${params.kingdom}.main.tsv.gz"), emit: main
  
    script:
    if (!params.taxonomy_aware) 
        """
        autometa-binning \
          --kmers $kmers \
          --coverages $coverage \
          --gc-content $gc_content \
          --markers $markers \
          --output-binning ${coverage.simpleName}.${params.kingdom}.binning.tsv.gz \
          --output-main ${coverage.simpleName}.${params.kingdom}.main.tsv.gz \
          --clustering-method ${params.clustering_method} \
          --completeness ${params.completeness} \
          --purity ${params.purity} \
          --cov-stddev-limit ${params.cov_stddev_limit} \
          --gc-stddev-limit ${params.gc_stddev_limit} \
          --starting-rank ${params.binning_starting_rank} \
          --domain ${params.kingdom}
        """
    else
        """
        autometa-binning \
          --kmers $kmers \
          --coverages $coverage \
          --gc-content $gc_content \
          --markers $markers \
          --output-binning ${coverage.simpleName}.${params.kingdom}.binning.tsv.gz \
          --output-main ${coverage.simpleName}.${params.kingdom}.main.tsv.gz \
          --clustering-method ${params.clustering_method} \
          --completeness ${params.completeness} \
          --purity ${params.purity} \
          --cov-stddev-limit ${params.cov_stddev_limit} \
          --gc-stddev-limit ${params.gc_stddev_limit} \
          --taxonomy $taxonomy \
          --starting-rank ${params.binning_starting_rank} \
          --domain ${params.kingdom}
        """
}
