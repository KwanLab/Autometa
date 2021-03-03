#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.length_cutoff = 3000
params.kmer_size = 4
params.cpus = 2
params.markers_database = "$HOME/Autometa/autometa/databases/markers"
params.kingdom = "bacteria"
params.metagenome = "</path/to/user/metagenome.fna>"
params.interim = "</path/to/store/user/interimediate/results>"
params.processed = "</path/to/store/user/final/results>"


process LENGTH_FILTER {
  tag "filtering metagenome ${metagenome.simpleName}"
  // container = 'placeholder for autometa image'
  publishDir params.interim, pattern: "${metagenome.simpleName}.filtered.fna", mode:'copy'

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.filtered.fna"

  """
  # autometa-length-filter [-h] [--cutoff <int>] [--force] [--verbose] [--stats] assembly out
  autometa-length-filter --cutoff ${params.length_cutoff} $metagenome ${metagenome.simpleName}.filtered.fna
  """
}

process KMERS {
  tag "counting kmers for ${metagenome.simpleName}"
  // container = 'placeholder for autometa image'
  cpus params.cpus
  publishDir params.processed, pattern: "*.kmers.*", mode:'copy'

  input:
    path metagenome

  output:
    path "*.kmers.tsv", emit: counts
    path "*.kmers.normalized.tsv", emit: normalized
    path "*.kmers.embedded.tsv", emit: embedded

  """
  #usage: autometa-kmers [-h] --fasta FASTA --kmers KMERS [--size SIZE] [--normalized NORMALIZED] [--norm-method {ilr,clr,am_clr}] [--do-pca] [--pca-dimensions PCA_DIMENSIONS] [--embedded EMBEDDED] [--embed-method {sksne,bhsne,umap}] [--embed-dimensions EMBED_DIMENSIONS] [--force] [--multiprocess] [--cpus CPUS] [--seed SEED]
  autometa-kmers \
    --fasta $metagenome \
    --kmers ${metagenome.simpleName}.kmers.tsv \
    --size ${params.kmer_size} \
    --normalized ${metagenome.simpleName}.kmers.normalized.tsv \
    --norm-method am_clr \
    --do-pca \
    --pca-dimensions 50 \
    --embedded ${metagenome.simpleName}.kmers.embedded.tsv \
    --embed-method bhsne \
    --embed-dimensions 2 \
    --multiprocess \
    --cpus ${task.cpus} \
    --seed 42
  """
}

process COVERAGE {
  tag "Calculating coverage for ${metagenome.simpleName}"
  // container = 'placeholder for autometa image'
  cpus params.cpus
  publishDir params.interim, pattern: "${metagenome.simpleName}.coverages.tsv", mode:'copy'

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.coverages.tsv"

  """
  # usage: autometa-coverage [-h] -f ASSEMBLY [-1 [FWD_READS [FWD_READS ...]]] [-2 [REV_READS [REV_READS ...]]] [-U [SE_READS [SE_READS ...]]] [--sam SAM] [--bam BAM] [--lengths LENGTHS] [--bed BED] [--cpus CPUS] [--from-spades] --out OUT
  autometa-coverage -f $metagenome --cpus ${task.cpus} --from-spades --out ${metagenome.simpleName}.coverages.tsv
  # Some other usage will have to be provided if the user needs to calculate coverage from reads. i.e.
  # autometa-coverage -f $metagenome --cpus ${task.cpus} -1 \$fwd_reads -2 \$rev_reads --out ${metagenome.simpleName}.coverages.tsv
  """
}

process MARKERS {
  tag "Finding markers for ${orfs.simpleName}"
  // container = 'placeholder for autometa image'
  cpus params.cpus
  publishDir params.interim, pattern: "${orfs.simpleName}.markers.tsv", mode:'copy'

  input:
    path orfs

  output:
    path "${orfs.simpleName}.markers.tsv"

  """
  # usage: autometa-markers [-h] [--orfs ORFS] [--kingdom {bacteria,archaea}] [--hmmscan HMMSCAN] [--out OUT] [--dbdir DBDIR] [--force] [--parallel] [--gnu-parallel] [--cpus CPUS] [--seed SEED]
  autometa-markers \
    --orfs $orfs \
    --out ${orfs.simpleName}.markers.tsv \
    --kingdom ${params.kingdom}
    --dbdir ${params.markers_database}
    --multiprocess \
    --cpus ${task.cpus} \
    --seed 42
  """
}

process ORFS {
  tag "Calling orfs for ${metagenome.simpleName}"
  // container = 'placeholder for autometa image'
  cpus params.cpus
  publishDir params.interim, pattern: '*.orfs.*', mode:'copy'

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.orfs.fna", emit: nucls
    path "${metagenome.simpleName}.orfs.faa", emit: prots

  """
  # usage: autometa-orfs [-h] [--force] [--cpus CPUS] [--parallel] assembly nucls_out prots_out
  autometa-orfs --cpus ${task.cpus} --parallel $metagenome ${metagenome.simpleName}.orfs.fna ${metagenome.simpleName}.orfs.faa
  """
}
