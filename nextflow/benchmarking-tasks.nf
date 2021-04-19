#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.community = "Community that was used for clustering or classification" //choices: 78, 156, 312, 625, 1250, 2500, 5000, 10000, etc.
params.processed = "Path to store final results"
params.ncbi_database = "Path to ncbi databases directory"


process BENCHMARK_CLUSTERING {
  tag "benchmarking clustering on ${community}"
  container = 'jason-c-kwan/autometa:dev'
  publishDir params.processed, pattern: "*.clustering_benchmarks.*"

  input:
    tuple path(binning), val(community)

  output:
    path "${community}.clustering_benchmarks.*.tsv.gz"

  script:
  """
  # First retrieve reference assignments for provided community
  autometa-datasets \
    --community ${community} \
    --file reference_assignments,binning \
    --output ${community}.reference.tsv.gz,${community}.binning.tsv.gz

  # Now benchmark inputs and previous gold-standard against community reference assignments
  autometa-benchmark \
    --benchmark clustering \
    --predictions $binning ${community}.binning.tsv.gz \
    --reference ${community}.reference.tsv.gz \
    --output-wide ${community}.clustering_benchmarks.wide.tsv.gz \
    --output-long ${community}.clustering_benchmarks.long.tsv.gz
  """
}

process BENCHMARK_CLASSIFICATION {
  tag "benchmarking classification on ${community}"
  container = 'jason-c-kwan/autometa:dev'
  containerOptions = "-v ${params.ncbi_database}:/ncbi:ro"
  publishDir params.processed, pattern: "*.classification_benchmarks.*"

  input:
    tuple path(taxonomy), val(community)

  output:
    path "${community}.classification_benchmarks.*.tsv.gz"

  script:
  """
  # First retrieve reference assignments for provided community
  autometa-datasets \
    --community ${community} \
    --file reference_assignments,taxonomy \
    --output ${community}.reference.tsv.gz,${community}.taxonomy.tsv.gz

  # Now benchmark inputs and previous gold-standard against community reference assignments
  autometa-benchmark \
    --benchmark classification \
    --predictions $taxonomy ${community}.taxonomy.tsv.gz \
    --reference ${community}.reference.tsv.gz \
    --output-wide ${community}.classification_benchmarks.wide.tsv.gz \
    --ncbi /ncbi
  """
}
