
params.fix_diamond_options      = [:]

// This workflow is used instead of DIAMOND_BLASTP directly to provide flexibility
// in possibly using other programs in the future

include { DIAMOND_BLASTP } from './../../modules/local/diamond_blastp.nf' addParams( options: params.diamond_blastp_options )

workflow SEARCH_TAXONOMY {
    take:
        orfs

    main:
        DIAMOND_BLASTP(orfs)

    emit:
        blastp_table = DIAMOND_BLASTP.out.blastp_table // BLAST fmt 6 table
}
