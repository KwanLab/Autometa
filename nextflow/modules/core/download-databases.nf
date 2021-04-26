#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.single_db_dir = "${projectDir}/large_db"

params.taxdump_dir = params.single_db_dir
params.accession2taxid_dir = params.single_db_dir
params.nr_dir = params.single_db_dir
params.bacteria_single_copy_dir = params.single_db_dir
params.bacteria_single_copy_cutoffs_dir = params.single_db_dir
params.archaea_single_copy_dir = params.single_db_dir
params.archaea_single_copy_cutoffs_dir = params.single_db_dir
 
include { SINGLE_RSYNC_DOWNLOAD as DOWNLOAD_TAXDUMP } from "../utilities/download.nf" addParams(rsync_storedir: params.taxdump_dir)
include { SINGLE_RSYNC_DOWNLOAD as DOWNLOAD_ACCESSION2TAXID_DIR } from "../utilities/download.nf" addParams(rsync_storedir: params.accession2taxid_dir)
include { SINGLE_RSYNC_DOWNLOAD as DOWNLOAD_NR } from "../utilities/download.nf" addParams(rsync_storedir: params.nr_dir)

include { SINGLE_WGET_DOWNLOAD as DOWNLOAD_BACTERIA_SINGLE_COPY } from "../utilities/download.nf" addParams(rsync_storedir: params.bacteria_single_copy_dir)
include { SINGLE_WGET_DOWNLOAD as DOWNLOAD_BACTERIA_SINGLE_COPY_CUTOFFS } from "../utilities/download.nf" addParams(rsync_storedir: params.bacteria_single_copy_cutoffs_dir)
include { SINGLE_WGET_DOWNLOAD as DOWNLOAD_ARCHAEA_SINGLE_COPY } from "../utilities/download.nf" addParams(rsync_storedir: params.archaea_single_copy_dir)
include { SINGLE_WGET_DOWNLOAD as DOWNLOAD_ARCHAEA_SINGLE_COPY_CUTOFFS } from "../utilities/download.nf" addParams(rsync_storedir: params.archaea_single_copy_cutoffs_dir)


workflow {
    DOWNLOAD_TAXDUMP(
        "rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        "taxdump.tar.gz"
    )
    DOWNLOAD_ACCESSION2TAXID_DIR(
        "rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz",
        "prot.accession2taxid.gz"
    )
    DOWNLOAD_NR(
        "rsync://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz",
        "nr.gz"
    )
    DOWNLOAD_BACTERIA_SINGLE_COPY(
        "https://raw.githubusercontent.com/KwanLab/Autometa/dev/autometa/databases/markers/bacteria.single_copy.hmm",
        "https://raw.githubusercontent.com/KwanLab/Autometa/dev/autometa/databases/markers/bacteria.single_copy.hmm.md5",
        "bacteria.single_copy.hmm"
    )
    DOWNLOAD_BACTERIA_SINGLE_COPY_CUTOFFS(
        "https://raw.githubusercontent.com/KwanLab/Autometa/dev/autometa/databases/markers/bacteria.single_copy.cutoffs",
        "https://raw.githubusercontent.com/KwanLab/Autometa/dev/autometa/databases/markers/bacteria.single_copy.cutoffs.md5",
        "bacteria.single_copy.cutoffs"
    )
    DOWNLOAD_ARCHAEA_SINGLE_COPY(
        "https://raw.githubusercontent.com/KwanLab/Autometa/dev/autometa/databases/markers/archaea.single_copy.hmm",
        "https://raw.githubusercontent.com/KwanLab/Autometa/dev/autometa/databases/markers/archaea.single_copy.hmm.md5",
        "archaea.single_copy.hmm"
    )
    DOWNLOAD_ARCHAEA_SINGLE_COPY_CUTOFFS(
        "https://raw.githubusercontent.com/KwanLab/Autometa/dev/autometa/databases/markers/archaea.single_copy.cutoffs",
        "https://raw.githubusercontent.com/KwanLab/Autometa/dev/autometa/databases/markers/archaea.single_copy.cutoffs.md5",
        "archaea.single_copy.cutoffs"
    )
}