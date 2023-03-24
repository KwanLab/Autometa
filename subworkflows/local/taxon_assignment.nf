
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

        // split channel based on the output fna files
        // add the fna file name as taxon to the meta map (e.g. "bacteria", "archaea")

        def retrieve_taxon_name  = {
            def filename = it.getName()
            // filename.split = [mock_data, archaea, fna, gz]
            filename.split("\\.")[1]
            }

        SPLIT_KINGDOMS
            .out
            .fasta
            .transpose()
            .multiMap{ meta, fasta ->
                   fasta:[meta + [taxon: retrieve_taxon_name(fasta)], fasta] // [[meta.id, meta.taxon], fasta]
                   taxon: retrieve_taxon_name(fasta)
                }
                .set{ taxon_split_fasta }

        // collect all the taxa from TAXON_ASSIGNMENT into a list (this should be variable)
        taxon_split_fasta
            .taxon
            .collect()
            .flatten()
            .set{ch_found_taxa_list}

        // TODO: not necessary but modifying autometa-taxonomy so that "taxonomy.tsv"
        // is split and named identically to the FASTA output would allow the logic here
        // to be identical to TAXON_ASSIGNMENT.out.taxon_split_fasta above
        // This creates a separate channel for each taxon but each with an identical SPLIT_KINGDOMS.out.contig_taxonomy_tsv
        SPLIT_KINGDOMS
            .out
            .contig_taxonomy_tsv
            .combine(ch_found_taxa_list)
            .map{ meta, tsv, taxon ->
                    [meta + [taxon: taxon], tsv]
                }
            .set{ contig_taxonomy_tsv }

    emit:
        contig_taxonomy_tsv = contig_taxonomy_tsv   // [[meta.id, meta.taxon], tsv]
        taxon_split_fasta   = taxon_split_fasta.fasta     // [[meta.id, meta.taxon], fasta]
        found_taxa_list     = ch_found_taxa_list    // [meta.taxon]
        orf_votes           = LCA.out.lca
        contig_votes        = MAJORITY_VOTE.out.votes
        taxdump_files       = PREPARE_TAXONOMY_DATABASES.out.taxdump_files
        versions            = ch_versions

}
