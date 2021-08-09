#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.community_type = "Community type that was used for clustering or classification" //choices: synthetic, simulated, all
params.community_size = "Community size that was used for clustering or classification" //choices: 78Mbp,156Mbp,312Mbp,625Mbp,1250Mbp,2500Mbp,5000Mbp,10000Mbp,all, etc.
params.processed = "Path to store final results"
params.ncbi_database = "Path to ncbi databases directory"


process BENCHMARK_CLUSTERING {
  tag "benchmarking clustering on ${community_type}: ${community_size}"
  container = 'jason-c-kwan/autometa:dev'
  publishDir params.processed, pattern: "*.clustering_benchmarks.*"

  input:
    tuple path(binning), val(community_type), val(community_size)

  output:
    path "${community_type}.${community_size}.clustering_benchmarks.*.tsv.gz"

  script:
  """
  # First retrieve reference assignments for provided community
  autometa-download-dataset \
    --community-type ${community_type} \
    --community-sizes ${community_size} \
    --file-names reference_assignments.tsv.gz,binning.tsv.gz \
    --dir-path .

  # Now benchmark inputs and previous gold-standard against community reference assignments
  autometa-benchmark \
    --benchmark clustering \
    --predictions $binning ${community_type}/${community_size}/binning.tsv.gz \
    --reference ${community_type}/${community_size}/reference_assignments.tsv.gz \
    --output-wide ${community_type}.${community_size}.clustering_benchmarks.wide.tsv.gz \
    --output-long ${community_type}.${community_size}.clustering_benchmarks.long.tsv.gz
  """
}

process BENCHMARK_CLASSIFICATION {
  tag "benchmarking classification on ${community_type}: ${community_size}"
  container = 'jason-c-kwan/autometa:dev'
  containerOptions = "-v ${params.ncbi_database}:/ncbi:ro"
  publishDir params.processed, pattern: "*.classification_benchmarks.*"

  input:
    tuple path(taxonomy), val(community_type), val(community_size)

  output:
    path "${community_type}.${community_size}.classification_benchmarks.*.tsv.gz"

  script:
  """
  # First retrieve reference assignments for provided community
  autometa-download-dataset \
    --community-type ${community_type} \
    --community-sizes ${community_size} \
    --file-names reference_assignments.tsv.gz,taxonomy.tsv.gz \
    --dir-path .

  # Now benchmark inputs and previous gold-standard against community reference assignments
  autometa-benchmark \
    --benchmark classification \
    --predictions $taxonomy ${community_type}/${community_size}/taxonomy.tsv.gz \
    --reference ${community_type}/${community_size}/reference_assignments.tsv.gz \
    --output-wide ${community_type}.${community_size}.classification_benchmarks.wide.tsv.gz \
    --ncbi /ncbi
  """
}
