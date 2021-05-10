#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process DIAMOND {
  label 'process_high'
  label 'process_long'
  
  tag "diamond blastp on ${orfs.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:ro"
  publishDir params.interim_dir, pattern: "${orfs.simpleName}.blastp.tsv"

  input:
    path orfs

  output:
    path "${orfs.simpleName}.blastp.tsv", emit: blastp_table

  """
  diamond blastp \
    --query ${orfs} \
    --db /ncbi/nr.dmnd \
    --evalue 1e-5 \
    --max-target-seqs 200 \
    --threads ${task.cpus} \
    --outfmt 6 \
    --out ${orfs.simpleName}.blastp.tsv \
    -b 10

  """
}