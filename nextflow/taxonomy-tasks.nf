#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.interim = "</path/to/store/user/interimediate/results>"
params.processed = "</path/to/store/user/final/results>"
params.ncbi_database = "$HOME/Autometa/autometa/databases/ncbi"
params.diamond_database = "$HOME/Autometa/autometa/databases/ncbi/nr.dmnd"
params.cpus = 2

process DIAMOND {
  tag "diamond blastp on ${orfs.simpleName}"
  // container = 'placeholder for autometa image'
  cpus params.cpus
  publishDir params.interim, pattern: "${orfs.simpleName}.blastp.tsv"

  input:
    path orfs

  output:
    path "${orfs.simpleName}.blastp.tsv"

  """
  diamond blastp \
    --query $orfs \
    --db ${params.diamond_database} \
    --evalue 1e-5 \
    --max-target-seqs 200 \
    --threads ${task.cpus} \
    --outfmt 6 \
    --out ${orfs.simpleName}.blastp.tsv
  """
}

process LCA {
  tag "Assigning LCA to ${orfs.simpleName}"
  // container = 'placeholder for autometa image'
  publishDir params.interim, pattern: "${orfs.simpleName}.lca.tsv"

  input:
    path orfs
    path blast

  output:
    path "${orfs.simpleName}.lca.tsv"

  """
  autometa-taxonomy-lca --blast $blast --dbdir ${params.ncbi_database} $orfs ${orfs.simpleName}.lca.tsv
  """
}

process MAJORITY_VOTE {
  tag "Performing taxon majority vote on ${orfs.simpleName}"
  // container = 'placeholder for autometa image'
  publishDir params.interim, pattern: "${orfs.simpleName}.votes.tsv"

  input:
    path orfs
    path lca

  output:
    path "${orfs.simpleName}.votes.tsv"

  """
  autometa-taxonomy-majority-vote --orfs $orfs --output ${orfs.simpleName}.votes.tsv --dbdir ${params.ncbi_database} --lca $lca
  """
}

process SPLIT_KINGDOMS {
  tag "Splitting votes into kingdoms for ${metagenome.simpleName}"
  // container = 'placeholder for autometa image'
  publishDir params.interim, pattern: "${metagenome.simpleName}.taxonomy.tsv"
  publishDir params.interim, pattern: '*.fna'

  input:
    path metagenome
    path votes

  output:
    path "${metagenome.simpleName}.taxonomy.tsv", emit: taxonomy
    path "${metagenome.simpleName}.bacteria.fna", emit: bacteria
    path "${metagenome.simpleName}.archaea.fna", emit: archaea

  """
  autometa-taxonomy \
    --input $votes \
    --output . \
    --prefix ${metagenome.simpleName} \
    --split-rank-and-write superkingdom \
    --assembly $metagenome \
    --ncbi ${params.ncbi_database}
  # Handling case where no archaea were recovered...
  if [[ ! -f ${metagenome.simpleName}.archaea.fna ]]
  then touch ${metagenome.simpleName}.archaea.fna
  fi
  """
}

// Autometa taxon assignment workflow
workflow TAXON_ASSIGNMENT {
    take:
      metagenome
      orfs

    main:
      DIAMOND(orfs)
      LCA(orfs, DIAMOND.out)
      MAJORITY_VOTE(orfs, LCA.out)
      SPLIT_KINGDOMS(metagenome, MAJORITY_VOTE.out)

    emit:
      taxonomy = SPLIT_KINGDOMS.out.taxonomy
      bacteria = SPLIT_KINGDOMS.out.bacteria
      archaea = SPLIT_KINGDOMS.out.archaea
      all_voting_results = LCA.out | mix(MAJORITY_VOTE.out) | collect
}
