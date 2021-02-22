#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.interim = "</path/to/store/user/interimediate/results>"
params.processed = "</path/to/store/user/final/results>"
params.ncbi_database = "$HOME/Autometa/autometa/databases/ncbi"
params.diamond_database = "$HOME/Autometa/autometa/databases/ncbi/nr.dmnd"
params.cpus = 2

process DIAMOND {
  tag "diamond blastp on ${orfs.simpleName}"
  container = 'placeholder for autometa image'
  cpus params.cpus
  publishDir params.interim, pattern: 'blastp.tsv', mode:'copy'

  input:
    path orfs

  output:
    path "blastp.tsv"

  """
  diamond blastp \
    --query $orfs \
    --db ${params.diamond_database} \
    --evalue 1e-5 \
    --max-target-seqs 200 \
    --threads ${task.cpus} \
    --outfmt 6 \
    --out blastp.tsv
  """
}

process LCA {
  tag "Assigning LCA to ${orfs.simpleName}"
  container = 'placeholder for autometa image'
  publishDir params.interim, pattern: 'lca.tsv', mode:'copy'

  input:
    path orfs
    path blast

  output:
    path "lca.tsv"

  """
  autometa-taxonomy-lca --blast $blast --dbdir ${params.ncbi_database} $orfs lca.tsv
  """
}

process MAJORITY_VOTE {
  tag "Performing taxon majority vote on ${metagenome.simpleName}"
  container = 'placeholder for autometa image'
  publishDir params.interim, pattern: 'votes.tsv', mode:'copy'

  input:
    path orfs
    path lca

  output:
    path "votes.tsv"

  """
  autometa-taxonomy-majority-vote --orfs $orfs --output votes.tsv --dbdir ${params.ncbi_database} --lca $lca
  """
}

process SPLIT_KINGDOMS {
  tag "Splitting votes into kingdoms for ${metagenome.simpleName}"
  container = 'placeholder for autometa image'
  publishDir params.processed, pattern: 'taxonomy.tsv', mode:'copy'
  publishDir params.processed, pattern: '*.fna', mode:'copy'

  input:
    path metagenome
    path orfs
    path votes

  output:
    path "taxonomy.tsv", emit: taxonomy
    path "Bacteria.fna", emit: bacteria
    path "Archaea.fna", emit: archaea
    // This may result in an error if there are not archaea present,
    // but I'm placing this here because we have written the pipeline to bin both bacteria and archaea

  """
  cp $votes votes.tsv.bkup
  autometa-taxonomy --split-rank-and-write superkingdom --assembly $metagenome --orfs $orfs --ncbi ${params.ncbi_database} $votes
  # script behavior should probably be cleaned up so we don't have to do these copies and renames.
  mv $votes taxonomy.tsv
  mv votes.tsv.bkup votes.tsv
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
      SPLIT_KINGDOMS(metagenome, orfs, MAJORITY_VOTE.out)

    emit:
      taxonomy = SPLIT_KINGDOMS.out.taxonomy
      bacteria = SPLIT_KINGDOMS.out.bacteria
      archaea = SPLIT_KINGDOMS.out.archaea
      all_voting_results = LCA.out | mix(MAJORITY_VOTE.out) | collect
}
