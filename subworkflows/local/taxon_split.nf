
include { PREPARE_NR_DB  } from './prepare_nr.nf'
include { PREPARE_TAXONOMY_DATABASES  } from './prepare_ncbi_taxinfo.nf'
include { LCA            } from './lca.nf'
include { MAJORITY_VOTE  } from './../../modules/local/majority_vote.nf'
include { SPLIT_KINGDOMS } from './../../modules/local/split_kingdoms.nf'
include { DIAMOND_BLASTP } from './../../modules/local/diamond_blastp.nf'




// Autometa taxon assignment workflow
workflow TAXON_SPLIT {
    take:
        contigs_and_orfs
        diamond_db_ch
        taxdump_ch
        prot_accession2taxid_ch
        dbtype_ch

    main:
        ch_versions = Channel.empty()

        contigs_and_orfs.multiMap { meta, fna_file, orfs_file ->
                fna: [meta, fna_file]
                orfs: [meta, orfs_file]
            }.set { result }

        DIAMOND_BLASTP (
            result.orfs,
            diamond_db_ch
        )
        ch_versions = ch_versions.mix(DIAMOND_BLASTP.out.versions)

        LCA (
            DIAMOND_BLASTP.out.diamond_results,
            taxdump_ch,
            prot_accession2taxid_ch,
            dbtype_ch
        )
        ch_versions = ch_versions.mix(LCA.out.versions)

        MAJORITY_VOTE (
            LCA.out.lca,
            taxdump_ch,
            dbtype_ch
        )
        ch_versions = ch_versions.mix(MAJORITY_VOTE.out.versions)

        result.fna
            .join(
                MAJORITY_VOTE.out.votes
            )
            .set{split_kingdoms_input}

        SPLIT_KINGDOMS (
            split_kingdoms_input,
            taxdump_ch,
            dbtype_ch
        )

        // Step 1: Generate combinations of meta and fna_file and flatten them correctly
        // handle if multiple fna files are present
        SPLIT_KINGDOMS.out.fna.map { meta, fna_file ->
            [[meta], fna_file].combinations()
        }.flatten().collate(2) // Creates pairs of [meta, fna_file]
        .set { tempch1 }

        // Step 2: Map each pair to set the taxon correctly for each meta-fna_file pair
        tempch1.map{  meta, fna_file ->
            // Set the taxon by extracting it from the fna_file name
            def new_meta = meta.clone()
            new_meta.taxon = fna_file.getName().tokenize('.')[-2]
            return [new_meta, fna_file] // Return a copy of meta to ensure independent taxon setting
        } .set { taxonomically_split_fna_ch }


        ch_versions = ch_versions.mix(SPLIT_KINGDOMS.out.versions)

    emit:
        taxonomy = SPLIT_KINGDOMS.out.taxonomy
        taxonomically_split_fna = taxonomically_split_fna_ch
        lca = LCA.out.lca
        votes = MAJORITY_VOTE.out.votes
        taxdump_files = taxdump_ch
        dbtype = dbtype_ch
        versions    = ch_versions

}
