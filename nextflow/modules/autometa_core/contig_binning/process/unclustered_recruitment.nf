#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process UNCLUSTERED_RECRUITMENT {
  tag "Performing Autometa unclustered recruitment"
  publishDir params.outdir, pattern: "${coverage.simpleName}.${params.kingdom}.recruitment.tsv.gz"

  input:
    path kmers
    path coverage
    path binning
    path markers
    val taxonomy

  output:
    path "${coverage.simpleName}.${params.kingdom}.recruitment.tsv.gz", emit: binning, optional: true 
    path "${coverage.simpleName}.${params.kingdom}.recruitment.main.tsv.gz", emit: main, optional: true 
  
  script:
  if (!params.taxonomy_aware) 
      """
      autometa-unclustered-recruitment \
        --classifier ${params.classification_method} \
        --kmer-dimensions ${params.classification_kmer_pca_dimensions} \
        --seed 42 \
        --kmers $kmers \
        --coverage $coverage \
        --binning $binning \
        --markers $markers \
        --output-binning ${coverage.simpleName}.${params.kingdom}.recruitment.tsv.gz \
        --output-main ${coverage.simpleName}.${params.kingdom}.recruitment.main.tsv.gz
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
         --output-binning ${coverage.simpleName}.${params.kingdom}.recruitment.tsv.gz \
         --output-main ${coverage.simpleName}.${params.kingdom}.recruitment.main.tsv.gz
       """
}
