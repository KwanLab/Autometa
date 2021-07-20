// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process MAJORITY_VOTE {
  label 'process_medium'

  tag "Performing taxon majority vote on ${meta.id}"
  publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

  conda (params.enable_conda ? "bioconda::autometa" : null)

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
       container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
  } else {
       container "jason-c-kwan/autometa:nfcore"
  }

  input:
    tuple val(meta), path(lca)
    path(ncbi_tax_dir)

  output:
    tuple val(meta), path("${meta.id}.votes.tsv"), emit: votes
    path  '*.version.txt'                        , emit: version

    script:
    // Add soft-links to original FastQs for consistent naming in pipeline
    def software = getSoftwareName(task.process)
  """
  autometa-taxonomy-majority-vote --lca ${lca} --output ${meta.id}.votes.tsv --dbdir "${ncbi_tax_dir}"

  echo "TODO" > autometa.version.txt
"""
}

