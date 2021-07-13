// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process ANALYZE_KMERS {
    tag "Counting kmers for ${meta.id}"
    label 'process_medium'
   
    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
   
    conda (params.enable_conda ? "autometa" : null)
 
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jason-c-kwan/autometa:nfcore"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
       tuple val(meta), path("${meta.id}.kmers.tsv"), emit: counts
       tuple val(meta), path("${meta.id}.kmers.normalized.tsv"), emit: normalized
       tuple val(meta), path("${meta.id}.kmers.embedded.tsv"), emit: embedded

    """
    autometa-kmers \
      --fasta $metagenome \
      --kmers ${meta.id}.kmers.tsv \
      --size ${params.kmer_size} \
      --norm-output ${meta.id}.kmers.normalized.tsv \
      --norm-method ${params.norm_method} \
      --pca-dimensions ${params.pca_dimensions} \
      --embedding-output ${meta.id}.kmers.embedded.tsv \
      --embedding-method ${params.embedding_method} \
      --embedding-dimensions ${params.embedding_dimensions} \
      --cpus ${task.cpus} \
      --seed 42
    """
}
