// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process UNCLUSTERED_RECRUITMENT {
    tag "Performing Autometa unclustered recruitment on ${meta.id}"
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
        tuple val(meta),path("?.kmers"), path("?.coverage"), path(binning), path("?.markers")
        val(taxonomy)

  output:
    tuple val(meta), path ("${meta.id}.${params.kingdom}.recruitment.tsv.gz"), emit: binning, optional: true 
    tuple val(meta), path ("${meta.id}.${params.kingdom}.recruitment.main.tsv.gz"), emit: main, optional: true 
  
  script:
  if (!params.taxonomy_aware) 
      """
      autometa-unclustered-recruitment \
        --classifier ${params.classification_method} \
        --kmer-dimensions ${params.classification_kmer_pca_dimensions} \
        --seed 42 \
        --kmers ?.kmers \
        --coverage ?.coverage \
        --binning $binning \
        --markers ?.markers \
        --output-binning ${meta.id}.${params.kingdom}.recruitment.tsv.gz \
        --output-main ${meta.id}.${params.kingdom}.recruitment.main.tsv.gz
      """
  else
       """
       autometa-unclustered-recruitment \
         --classifier ${params.classification_method} \
         --kmer-dimensions ${params.classification_kmer_pca_dimensions} \
         --seed 42 \
         --taxonomy $taxonomy \
         --kmers $kmers \
         --coverage $coverage \
         --binning $binning \
         --markers $markers \
         --output-binning ${meta.id}.${params.kingdom}.recruitment.tsv.gz \
         --output-main ${meta.id}.${params.kingdom}.recruitment.main.tsv.gz
       """
}
