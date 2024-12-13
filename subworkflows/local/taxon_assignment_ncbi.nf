
include { PREPARE_NR_DB                 } from './prepare_nr.nf'
include { PREPARE_TAXONOMY_DATABASES    } from './prepare_ncbi_taxinfo.nf'
include { TAXON_SPLIT                   } from './taxon_split.nf'

// Autometa taxon assignment workflow
workflow NCBI_TAXON_ASSIGNMENT {
    take:
        filtered_metagenome_fasta
        orfs

    main:
        ch_versions = Channel.empty()
        dbtype_ch = channel.value( 'ncbi')

        PREPARE_TAXONOMY_DATABASES()
        ch_versions = ch_versions.mix(PREPARE_TAXONOMY_DATABASES.out.versions)

        PREPARE_NR_DB()
        ch_versions = ch_versions.mix(PREPARE_NR_DB.out.versions)

        contigs_and_orfs = filtered_metagenome_fasta.join(orfs)

        TAXON_SPLIT(
            contigs_and_orfs,
            PREPARE_NR_DB.out.diamond_db,
            PREPARE_TAXONOMY_DATABASES.out.taxdump_files,
            PREPARE_TAXONOMY_DATABASES.out.prot_accession2taxid,
            dbtype_ch
        )

        ch_versions = ch_versions.mix(TAXON_SPLIT.out.versions)

    emit:
        taxonomy                = TAXON_SPLIT.out.taxonomy
        taxonomically_split_fna = TAXON_SPLIT.out.taxonomically_split_fna
        lca                     = TAXON_SPLIT.out.lca
        votes                   = TAXON_SPLIT.out.votes
        taxdump_files           = PREPARE_TAXONOMY_DATABASES.out.taxdump_files
        dbtype                  = dbtype_ch
        versions                = ch_versions

}

