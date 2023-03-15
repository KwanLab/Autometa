
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
        ch_versions = Channel.empty()

        PREPARE_TAXONOMY_DATABASES()
        ch_versions = ch_versions.mix(PREPARE_TAXONOMY_DATABASES.out.versions)

        PREPARE_NR_DB()
        ch_versions = ch_versions.mix(PREPARE_NR_DB.out.versions)

        DIAMOND_BLASTP (
            merged_prodigal,
            PREPARE_NR_DB.out.diamond_db
        )
        ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

        LCA (
            DIAMOND_BLASTP.out.diamond_results,
            PREPARE_TAXONOMY_DATABASES.out.taxdump_files,
            PREPARE_TAXONOMY_DATABASES.out.prot_accession2taxid
        )
        ch_versions = ch_versions.mix(LCA.out.versions)

        MAJORITY_VOTE (
            LCA.out.lca,
            PREPARE_TAXONOMY_DATABASES.out.taxdump_files
        )
        ch_versions = ch_versions.mix(MAJORITY_VOTE.out.versions)

        filtered_metagenome_fasta
            .join(
                MAJORITY_VOTE.out.votes
            )
            .set{split_kingdoms_input}

        SPLIT_KINGDOMS (
            split_kingdoms_input,
            PREPARE_TAXONOMY_DATABASES.out.taxdump_files
        )
        ch_versions = ch_versions.mix(SPLIT_KINGDOMS.out.versions)

    emit:
        taxonomy            = SPLIT_KINGDOMS.out.taxonomy
        taxon_split_fasta   = SPLIT_KINGDOMS.out.fasta
        orf_votes           = LCA.out.lca
        contig_votes        = MAJORITY_VOTE.out.votes
        taxdump_files       = PREPARE_TAXONOMY_DATABASES.out.taxdump_files
        versions            = ch_versions

}
