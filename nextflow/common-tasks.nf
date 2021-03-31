#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.length_cutoff = 3000
params.kmer_size = 5
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
  publishDir params.interim, pattern: "${metagenome.simpleName}.*"

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.filtered.fna", emit: fasta
    path "${metagenome.simpleName}.stats.tsv", emit: stats
    path "${metagenome.simpleName}.gc_content.tsv", emit: gc_content

  """
  autometa-length-filter \
    --assembly $metagenome \
    --cutoff ${params.length_cutoff} \
    --output-fasta ${metagenome.simpleName}.filtered.fna \
    --output-stats ${metagenome.simpleName}.stats.tsv \
    --output-gc-content ${metagenome.simpleName}.gc_content.tsv
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
    --cpus ${task.cpus} \
    --seed 42
  """
}

process KMER_COVERAGE {
  tag "Calculating k-mer coverage for ${metagenome.simpleName}"
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
  // copying orfs via stageInMode is required to run hmmscan (does not handle symlinks)
  stageInMode 'copy'
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

process ALIGN_READS {
  tag "Aligning reads to ${metagenome.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  cpus params.cpus

  input:
    path metagenome
    path fwd_reads
    path rev_reads
    path se_reads

  output:
    path "${metagenome.simpleName}.sam"

  """
  bowtie2-build \
    --threads ${task.cpus}
    ${metagenome} \
    ${metagenome.simpleName}.db

  bowtie2 \
    -x ${metagenome.simpleName}.db \
    -q \
    --phred33 \
    --very-sensitive \
    --no-unal \
    -p ${task.cpus} \
    -S ${metagenome.simpleName}.sam \
    -1 $fwd_reads \
    -2 $rev_reads \
    -U $se_reads
  """
}

process SORT_READS {
  tag "Sorting reads to ${sam.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  cpus params.cpus

  input:
    path sam

  output:
    path "${sam.simpleName}.bam"

  """
  samtools view -@${task.cpus} -bS ${sam} \
    | samtools sort -@${task.cpus} -o ${sam.simpleName}.bam
  """
}

process LENGTH_TABLE {
  tag "length table for ${metagenome.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  cpus params.cpus

  input:
    path metagenome

  output:
    path "${metagenome.simpleName}.lengths.tsv"

  """
  #!/usr/bin/env python

  from Bio import SeqIO
  import pandas as pd

  seqs = {record.id: len(record) for record in SeqIO.parse($metagenome, "fasta")}
  lengths = pd.Series(seqs, name="length")
  lengths.index.name = "contig"
  lengths.to_csv(${metagenome.simpleName}.lengths.tsv, sep="\t", index=True, header=True)
  """
}

process GENOMECOV {
  tag "Computing genome coverage for ${bam.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  cpus params.cpus

  input:
    path bam
    path lengths

  output:
    path "${bam.simpleName}.bed.tsv", emit: bed
    path "${bam.simpleName}.coverage.tsv", emit: coverage

  """
  bedtools genomecov -ibam $bam -g $lengths > ${bam.simpleName}.bed.tsv
  autometa-parse-bed \
    --ibam $bam \
    --lengths $lengths \
    --bed ${bam.simpleName}.bed.tsv \
    --output coverage.tsv
  """
}

workflow READ_COVERAGE {
  take:
    metagenome
    fwd_reads
    rev_reads
    se_reads

  main:
    LENGTH_TABLE(metagenome)
    ALIGN_READS(metagenome, fwd_reads, rev_reads, se_reads)
    SORT_READS(ALIGN_READS.out)
    GENOMECOV(SORT_READS.out, LENGTH_TABLE.out)

  emit:
    bed = GENOMECOV.out.bed
    coverage = GENOMECOV.out.coverage
}
