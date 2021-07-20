params.lca_options = [:]
params.majority_vote_options = [:]
params.split_kingdoms_options = [:]

include { LCA } from './../../modules/local/lca.nf' addParams( options: params.lca_options )
include { MAJORITY_VOTE } from './../../modules/local/majority_vote.nf' addParams( options: params.majority_vote_options )
include { SPLIT_KINGDOMS } from './../../modules/local/split_kingdoms.nf' addParams( options: params.split_kingdoms_options )


// Autometa taxon assignment workflow
workflow TAXON_ASSIGNMENT {
    take:
        metagenome
        blastp_table
        ncbi_tax_dir

    main:
        LCA(blastp_table, ncbi_tax_dir) // output '${blast.simpleName}.lca.tsv'
        MAJORITY_VOTE(LCA.out.lca, ncbi_tax_dir) //output ${lca.simpleName}.votes.tsv

        metagenome
        .join(
            MAJORITY_VOTE.out.votes
        )
        .set{split_kingdoms_input}


        SPLIT_KINGDOMS(split_kingdoms_input, ncbi_tax_dir)
        // output "${assembly.simpleName}.taxonomy.tsv" "${assembly.simpleName}.bacteria.fna" , "${assembly.simpleName}.archaea.fna"

    emit:
        taxonomy = SPLIT_KINGDOMS.out.taxonomy
        bacteria = SPLIT_KINGDOMS.out.bacteria
        archaea = SPLIT_KINGDOMS.out.archaea
        orf_votes = LCA.out.lca
        contig_votes = MAJORITY_VOTE.out.votes
}


