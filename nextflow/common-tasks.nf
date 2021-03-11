#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.length_cutoff = 3000
params.kmer_size = 4
params.kmer_norm_method = "am_clr" // choices: "am_clr", "clr", "ilr"
params.kmer_pca_dimensions = 50
params.kmer_embed_method = "bhsne" // choices: "sksne", "bhsne", "umap"
params.kmer_embed_dimensions = 2
params.cpus = 2
params.kingdom = "bacteria"
params.metagenome = "</path/to/user/metagenome.fna>"
params.interim = "</path/to/store/user/interimediate/results>"
params.processed = "</path/to/store/user/final/results>"


process LENGTH_FILTER {
  tag "filtering metagenome ${metagenome.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  publishDir params.interim, pattern: "${metagenome.simpleName}.filtered.fna"

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.filtered.fna"

  """
  autometa-length-filter --cutoff ${params.length_cutoff} $metagenome ${metagenome.simpleName}.filtered.fna
  """
}

process KMERS {
  tag "counting kmers for ${metagenome.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  cpus params.cpus
  publishDir params.interim, pattern: "*.kmers.*"

  input:
    path metagenome

  output:
    path "*.kmers.tsv", emit: counts
    path "*.kmers.normalized.tsv", emit: normalized
    path "*.kmers.embedded.tsv", emit: embedded

  """
  autometa-kmers \
    --fasta $metagenome \
    --kmers ${metagenome.simpleName}.kmers.tsv \
    --size ${params.kmer_size} \
    --normalized ${metagenome.simpleName}.kmers.normalized.tsv \
    --norm-method ${params.kmer_norm_method} \
    --do-pca \
    --pca-dimensions ${params.kmer_pca_dimensions} \
    --embedded ${metagenome.simpleName}.kmers.embedded.tsv \
    --embed-method bhsne \
    --embed-dimensions ${params.kmer_embed_dimensions} \
    --multiprocess \
    --cpus ${task.cpus} \
    --seed 42
  """
}

process COVERAGE {
  tag "Calculating coverage for ${metagenome.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  cpus params.cpus
  publishDir params.interim, pattern: "${metagenome.simpleName}.coverages.tsv"

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.coverages.tsv"

  """
  autometa-coverage \
    --assembly $metagenome \
    --cpus ${task.cpus} \
    --from-spades \
    --out ${metagenome.simpleName}.coverages.tsv
  """
}

process MARKERS {
  tag "Finding markers for ${orfs.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  cpus params.cpus
  publishDir params.interim, pattern: "${orfs.simpleName}.markers.tsv"
  publishDir params.interim, pattern: "${orfs.simpleName}.hmmscan.tsv"

  input:
    path orfs

  output:
    path "${orfs.simpleName}.markers.tsv"

  """
  autometa-markers \
    --orfs $orfs \
    --hmmscan ${orfs.simpleName}.hmmscan.tsv \
    --out ${orfs.simpleName}.markers.tsv \
    --kingdom ${params.kingdom} \
    --dbdir /Autometa/autometa/databases/markers \
    --parallel \
    --cpus ${task.cpus} \
    --seed 42
  """
}

process ORFS {
  tag "Calling orfs for ${metagenome.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  // Hardcoding cpus here b/c prodigal is limited to only using single core
  cpus 1
  publishDir params.interim, pattern: "${metagenome.simpleName}.orfs.f*"

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.orfs.fna", emit: nucls
    path "${metagenome.simpleName}.orfs.faa", emit: prots

  """
  prodigal \
    -i $metagenome \
    -d ${metagenome.simpleName}.orfs.fna \
    -a ${metagenome.simpleName}.orfs.faa \
    -p meta \
    -q \
    -m
  """
}
