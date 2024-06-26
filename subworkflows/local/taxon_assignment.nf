params.prepare_lca_options  = [:]
params.reduce_lca_options    = [:]
params.majority_vote_options = [:]
params.split_kingdoms_options = [:]
params.nr_dmnd_dir = [:]
params.taxdump_tar_gz_dir = [:]
params.prot_accession2taxid_gz_dir = [:]
params.diamond_blastp_options = [:]

params.debug = [:]
params.diamond_makedb_options = [:]
params.large_downloads_permission = [:]


include { PREPARE_NR_DB  } from './prepare_nr.nf'                         addParams( debug: params.debug, diamond_makedb_options: params.diamond_makedb_options, nr_dmnd_dir: params.nr_dmnd_dir  )
include { PREPARE_TAXONOMY_DATABASES  } from './prepare_ncbi_taxinfo.nf'  addParams( debug: params.debug, taxdump_tar_gz_dir: params.taxdump_tar_gz_dir, prot_accession2taxid_gz_dir: params.prot_accession2taxid_gz_dir  )
include { LCA            } from './lca.nf'                                addParams( prepare_lca_options: params.prepare_lca_options, reduce_lca_options: params.reduce_lca_options )
include { MAJORITY_VOTE  } from './../../modules/local/majority_vote.nf'  addParams( options: params.majority_vote_options       )
include { SPLIT_KINGDOMS } from './../../modules/local/split_kingdoms.nf' addParams( options: params.split_kingdoms_options      )
include { DIAMOND_BLASTP } from './../../modules/local/diamond_blastp.nf' addParams( options: params.diamond_blastp_options      )


// Autometa taxon assignment workflow
workflow TAXON_ASSIGNMENT {
    take:
        metagenome
        merged_prodigal

    main:
        // check if user has given permission for large downloads
        if (params.large_downloads_permission) {
            // Download and prep necessary databases
            PREPARE_NR_DB()
            PREPARE_NR_DB.out.diamond_db
                .set{diamond_db}
            PREPARE_TAXONOMY_DATABASES()
            PREPARE_TAXONOMY_DATABASES.out.taxdump
                .set{ncbi_taxdump}
            PREPARE_TAXONOMY_DATABASES.out.prot_accession2taxid
                .set{prot_accession2taxid}
        } else {
            // check for nr.dmnd, if not found, check for nr.gz
            // if nr.gz exists, create nr.dmnd
            // if nr.gz also doesn't exist, stop the pipeline
            if (!file("${params.nr_dmnd_dir}/nr.dmnd").exists()) {
                if (file("${params.nr_dmnd_dir}/nr.gz").exists()) {
                    PREPARE_NR_DB()
                    PREPARE_NR_DB.out.diamond_db
                        .set{diamond_db}
                } else {
                    throw new Exception("Neither nr.dmnd or nr.gz was found")
                }
            } else {
                diamond_db = file("${params.nr_dmnd_dir}/nr.dmnd", checkIfExists: true)
            }
        }

        DIAMOND_BLASTP (
            merged_prodigal,
            diamond_db
        )

        ncbi_tax_dir = file(params.taxdump_tar_gz_dir)

        LCA (
            DIAMOND_BLASTP.out.diamond_results,
            ncbi_tax_dir
        ) // output '${blast.simpleName}.lca.tsv'

        MAJORITY_VOTE (
            LCA.out.lca,
            ncbi_tax_dir
        ) //output ${lca.simpleName}.votes.tsv

        metagenome
            .join(
                MAJORITY_VOTE.out.votes
            )
            .set{split_kingdoms_input}

        SPLIT_KINGDOMS (
            split_kingdoms_input,
            ncbi_tax_dir
        )

    emit:
        taxonomy = SPLIT_KINGDOMS.out.taxonomy
        bacteria = SPLIT_KINGDOMS.out.bacteria
        archaea = SPLIT_KINGDOMS.out.archaea
        orf_votes = LCA.out.lca
        contig_votes = MAJORITY_VOTE.out.votes
}
