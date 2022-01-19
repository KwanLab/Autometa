#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.prepare_lca_options  = [:]
params.reduce_lca_options    = [:]

include { PREPARE_LCA    } from './../../modules/local/prepare_lca.nf' addParams( options: params.prepare_lca_options )
include { REDUCE_LCA     } from './../../modules/local/reduce_lca.nf'  addParams( options: params.reduce_lca_options )


workflow LCA {

    take:
        blastp_results
        blastp_dbdir

    main:
        PREPARE_LCA(
            blastp_dbdir
        )
        REDUCE_LCA(
            blastp_results,
            blastp_dbdir,
            PREPARE_LCA.out.cache
        )

    emit:
        lca = REDUCE_LCA.out.lca
        error_taxid = REDUCE_LCA.out.error_taxids
        sseqid_to_taxids = REDUCE_LCA.out.sseqid_to_taxids
        cache  = PREPARE_LCA.out.cache
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
