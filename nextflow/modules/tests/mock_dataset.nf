#!/usr/bin/env nextflow
nextflow.enable.dsl = 2



include { GET_ASSEMBLY_SUMMARY } from './process/download_mock_data.nf'
include { GET_FTP_DIRS } from './process/download_mock_data.nf'
include { DOWNLOAD_MOCK_DATA } from './process/download_mock_data.nf'
include { WRITE_FILE } from './process/download_mock_data.nf'


include { CREATE_DIAMOND_DB } from './process/create_blast_db.nf'


assemblies = Channel.fromList(
    [
        "GCF_008124965.1",
        "GCF_016882605.1",
        "GCF_900099615.1",
        "GCF_014654435.1",        
    ]
)

workflow DL {

    main:

    GET_ASSEMBLY_SUMMARY()
    GET_FTP_DIRS(
        GET_ASSEMBLY_SUMMARY.out,
        assemblies.flatten()
        )
    DOWNLOAD_MOCK_DATA(GET_FTP_DIRS.out)
    CREATE_DIAMOND_DB( DOWNLOAD_MOCK_DATA.out.protein.splitFasta(by:1).collectFile() )
    WRITE_FILE(
        DOWNLOAD_MOCK_DATA.out.nucleotide.splitFasta(by:1).collectFile(),
        "mock_metagenome.fna"
    )
}

workflow { DL()}
