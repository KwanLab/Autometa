#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.interim_dir = "</path/to/store/user/interimediate/results>"
params.outdir = "</path/to/store/user/final/results>"
params.ncbi_database = "$HOME/Autometa/autometa/databases/ncbi"
params.cpus = 2

process DIAMOND {
  tag "diamond blastp on ${orfs.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:ro"
  cpus params.cpus
  publishDir params.interim_dir, pattern: "${orfs.simpleName}.blastp.tsv"

  input:
    path orfs

  output:
    path "${orfs.simpleName}.blastp.tsv"

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

process LCA {
  tag "Assigning LCA to ${blast.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:rw"
  publishDir params.interim_dir, pattern: "${blast.simpleName}.lca.tsv"

  input:
    path blast

  output:
    path "${blast.simpleName}.lca.tsv"

  """
  autometa-taxonomy-lca --blast ${blast} --dbdir /ncbi --output ${blast.simpleName}.lca.tsv
  """
}

process MAJORITY_VOTE {
  tag "Performing taxon majority vote on ${lca.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:rw"
  publishDir params.interim_dir, pattern: "${lca.simpleName}.votes.tsv"

  input:
    path lca

  output:
    path "${lca.simpleName}.votes.tsv"

  """
  autometa-taxonomy-majority-vote --lca ${lca} --output ${lca.simpleName}.votes.tsv --dbdir /ncbi
  """
}

process SPLIT_KINGDOMS {
  tag "Splitting votes into kingdoms for ${assembly.simpleName}"
  containerOptions = "-v ${params.single_db_dir}:/ncbi:rw"
  publishDir params.interim_dir, pattern: "${assembly.simpleName}.taxonomy.tsv"
  publishDir params.interim_dir, pattern: '*.{bacteria,archaea}.fna'

  input:
    path votes
    path assembly

  output:
    path "${assembly.simpleName}.taxonomy.tsv", emit: taxonomy
    path "${assembly.simpleName}.bacteria.fna", emit: bacteria
    path "${assembly.simpleName}.archaea.fna", emit: archaea

  """
  autometa-taxonomy \
    --votes ${votes} \
    --output . \
    --prefix ${assembly.simpleName} \
    --split-rank-and-write superkingdom \
    --assembly ${assembly} \
    --ncbi /ncbi
  # Handling case where no archaea were recovered...
  if [[ ! -f ${assembly.simpleName}.archaea.fna ]]
  then touch ${assembly.simpleName}.archaea.fna
  fi
  """
}

// Autometa taxon assignment workflow
workflow TAXON_ASSIGNMENT {
    take:
      assembly
      orfs

    main:
      DIAMOND(orfs)
      LCA(DIAMOND.out)
      MAJORITY_VOTE(LCA.out)
      SPLIT_KINGDOMS(MAJORITY_VOTE.out, assembly)

    emit:
      taxonomy = SPLIT_KINGDOMS.out.taxonomy
      bacteria = SPLIT_KINGDOMS.out.bacteria
      archaea = SPLIT_KINGDOMS.out.archaea
      orf_votes = LCA.out
      contig_votes = MAJORITY_VOTE.out
}
