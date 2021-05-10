#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process BINNING_SUMMARY {
  tag "Binning summary for ${binning_main.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:ro"

  input:
    path binning_main
    path markers
    path metagenome
    val binning_column

  output:
    path 'metabin_stats.tsv', emit: stats
    path 'metabin_taxonomy.tsv', emit: taxonomies
    path 'metabins', emit: metabins

  script:
  """
  autometa-binning-summary \
    --ncbi /ncbi \
    --binning-main $binning_main \
    --markers $markers \
    --metagenome $metagenome \
    --binning-column $binning_column \
    --output-stats metabin_stats.tsv \
    --output-taxonomy metabin_taxonomy.tsv \
    --output-metabins metabins
  """
}
