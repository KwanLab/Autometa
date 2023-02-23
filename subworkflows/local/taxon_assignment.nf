
include { PREPARE_NR_DB  } from './prepare_nr.nf'
include { PREPARE_TAXONOMY_DATABASES  } from './prepare_ncbi_taxinfo.nf'
include { LCA            } from './lca.nf'
include { MAJORITY_VOTE  } from './../../modules/local/majority_vote.nf'
include { SPLIT_KINGDOMS } from './../../modules/local/split_kingdoms.nf'
include { DIAMOND_BLASTP } from './../../modules/local/diamond_blastp.nf'


// Autometa taxon assignment workflow
workflow TAXON_ASSIGNMENT {
    take:
        filtered_metagenome_fasta
        merged_prodigal

    main:

        PREPARE_TAXONOMY_DATABASES()
        PREPARE_NR_DB()

        DIAMOND_BLASTP (
            merged_prodigal,
            PREPARE_NR_DB.out.diamond_db
        )
        LCA (
            DIAMOND_BLASTP.out.diamond_results,
            PREPARE_TAXONOMY_DATABASES.out.taxdump_files,
            PREPARE_TAXONOMY_DATABASES.out.prot_accession2taxid
        ) // output '${blast.simpleName}.lca.tsv'

        MAJORITY_VOTE (
            LCA.out.lca,
            PREPARE_TAXONOMY_DATABASES.out.taxdump_files
        ) //output ${lca.simpleName}.votes.tsv

        filtered_metagenome_fasta
            .join(
                MAJORITY_VOTE.out.votes
            )
            .set{split_kingdoms_input}

        SPLIT_KINGDOMS (
            split_kingdoms_input,
            PREPARE_TAXONOMY_DATABASES.out.taxdump_files
        )

    emit:
        taxonomy = SPLIT_KINGDOMS.out.taxonomy
        bacteria = SPLIT_KINGDOMS.out.bacteria
        archaea = SPLIT_KINGDOMS.out.archaea
        orf_votes = LCA.out.lca
        contig_votes = MAJORITY_VOTE.out.votes
        taxdump_files = PREPARE_TAXONOMY_DATABASES.out.taxdump_files
}
