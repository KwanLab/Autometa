#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { PREPARE_LCA as PREP_DBS } from './../../modules/local/prepare_lca.nf'
include { REDUCE_LCA as REDUCE    } from './../../modules/local/reduce_lca.nf'


workflow LCA {

    take:
        blastp_results
        taxdump_files
        prot_accession2taxid

    main:

        PREP_DBS(
            taxdump_files
        )
        REDUCE(
            blastp_results,
            taxdump_files,
            PREP_DBS.out.cache,
            prot_accession2taxid
        )

    emit:
        lca = REDUCE.out.lca
        error_taxid = REDUCE.out.error_taxids
        sseqid_to_taxids = REDUCE.out.sseqid_to_taxids
        cache  = PREP_DBS.out.cache
}

/*
---------------: Test Command :-----------------
nextflow run -resume $HOME/Autometa/subworkflows/local/lca.nf \\
    --publish_dir_mode copy \\
    --input '/path/to/*.blastp.tsv' \\
    --dbdir '/path/to/ncbi/databases/directory'
*/

workflow {
    blastp_results_ch = Channel
            .fromPath(params.input)
            .map { row ->
                    def meta = [:]
                    meta.id = row.simpleName
                    return [ meta, row ]
                }

    dbdir = Channel.value(params.dbdir)

    LCA(
        blastp_results_ch,
        dbdir
    )
}
