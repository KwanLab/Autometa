#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.interim = "</path/to/store/user/interimediate/results>"
params.processed = "</path/to/store/user/final/results>"
params.ncbi_database = "$HOME/Autometa/autometa/databases/ncbi"
params.cpus = 2

process DIAMOND {
  tag "diamond blastp on ${orfs.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  containerOptions = "-v ${params.ncbi_database}:/ncbi:ro"
  cpus params.cpus
  publishDir params.interim, pattern: "${orfs.simpleName}.blastp.tsv"

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
    --out ${orfs.simpleName}.blastp.tsv
  """
}

process LCA {
  tag "Assigning LCA to ${orfs.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  containerOptions = "-v ${params.ncbi_database}:/ncbi:rw"
  publishDir params.interim, pattern: "${orfs.simpleName}.lca.tsv"

  input:
    path orfs
    path blast

  output:
    path "${orfs.simpleName}.lca.tsv"

  """
  autometa-taxonomy-lca --blast $blast --dbdir /ncbi ${orfs} ${orfs.simpleName}.lca.tsv
  """
}

process MAJORITY_VOTE {
  tag "Performing taxon majority vote on ${orfs.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  containerOptions = "-v ${params.ncbi_database}:/ncbi:rw"
  publishDir params.interim, pattern: "${orfs.simpleName}.votes.tsv"

  input:
    path orfs
    path lca

  output:
    path "${orfs.simpleName}.votes.tsv"

  """
  autometa-taxonomy-majority-vote --orfs ${orfs} --output ${orfs.simpleName}.votes.tsv --dbdir /ncbi --lca $lca
  """
}

process SPLIT_KINGDOMS {
  tag "Splitting votes into kingdoms for ${assembly.simpleName}"
  container = 'jason-c-kwan/autometa:dev'
  containerOptions = "-v ${params.ncbi_database}:/ncbi:rw"
  publishDir params.interim, pattern: "${assembly.simpleName}.taxonomy.tsv"
  publishDir params.interim, pattern: '*.{bacteria,archaea}.fna'

  input:
    path votes
    path assembly

  output:
    path "${assembly.simpleName}.taxonomy.tsv", emit: taxonomy
    path "${assembly.simpleName}.bacteria.fna", emit: bacteria
    path "${assembly.simpleName}.archaea.fna", emit: archaea

  """
  autometa-taxonomy \
    --input ${votes} \
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
      LCA(orfs, DIAMOND.out)
      MAJORITY_VOTE(orfs, LCA.out)
      SPLIT_KINGDOMS(MAJORITY_VOTE.out, assembly)

    emit:
      taxonomy = SPLIT_KINGDOMS.out.taxonomy
      bacteria = SPLIT_KINGDOMS.out.bacteria
      archaea = SPLIT_KINGDOMS.out.archaea
      orf_votes = LCA.out
      contig_votes = MAJORITY_VOTE.out
}
