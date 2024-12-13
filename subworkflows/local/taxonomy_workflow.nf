
include { NCBI_TAXON_ASSIGNMENT as NCBI  } from './taxon_assignment_ncbi.nf'
include { GTDB_TAXON_ASSIGNMENT as GTDB_REFINEMENT  } from './taxon_assignment_gtdb.nf'

// Autometa taxon assignment workflow
workflow TAXONOMY_WORKFLOW {
    take:
        filtered_metagenome_fasta
        merged_prodigal

    main:
        ch_versions = Channel.empty()

        if (params.taxonomy_aware) {
            NCBI(
                filtered_metagenome_fasta,
                merged_prodigal
            )
            ch_versions = ch_versions.mix(NCBI.out.versions)

            if (params.use_gtdb) {
                GTDB_REFINEMENT(
                    NCBI.out.taxonomically_split_fna,
                    merged_prodigal
                )
                ch_versions = ch_versions.mix(GTDB_REFINEMENT.out.versions)
                taxonomy = GTDB_REFINEMENT.out.taxonomy
                taxonomically_split_fna_ch = GTDB_REFINEMENT.out.taxonomically_split_fna
                orf_votes = GTDB_REFINEMENT.out.lca
                contig_votes = GTDB_REFINEMENT.out.votes
                taxdump_files = GTDB_REFINEMENT.out.taxdump_files
                dbtype = GTDB_REFINEMENT.out.dbtype

            } else {
                taxonomy = NCBI.out.taxonomy
                taxonomically_split_fna_ch = NCBI.out.taxonomically_split_fna
                orf_votes = NCBI.out.lca
                contig_votes = NCBI.out.votes
                taxdump_files = NCBI.out.taxdump_files
                dbtype = NCBI.out.dbtype
            }
     }

    emit:
        taxonomy = taxonomy
        taxonomically_split_fna = taxonomically_split_fna_ch
        orf_votes = orf_votes
        contig_votes = contig_votes
        taxdump_files = taxdump_files
        dbtype = dbtype
        versions    = ch_versions
}

