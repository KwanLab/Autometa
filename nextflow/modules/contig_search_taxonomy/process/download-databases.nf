#!/usr/bin/env nextflow
nextflow.enable.dsl=2


params.single_db_dir = null

// Allows flexibility for storing databases
// Defaults to everything saving into one directory, as in Autometa 1.0
params.taxdump_dir = params.single_db_dir
params.accession2taxid_dir = params.single_db_dir
params.nr_dir = params.single_db_dir
params.bacteria_single_copy_dir = params.single_db_dir
params.bacteria_single_copy_cutoffs_dir = params.single_db_dir
params.archaea_single_copy_dir = params.single_db_dir
params.archaea_single_copy_cutoffs_dir = params.single_db_dir
 
include { SINGLE_RSYNC_DOWNLOAD as DOWNLOAD_TAXDUMP } from "../utilities/download.nf" addParams(rsync_storedir: params.taxdump_dir)
include { SINGLE_RSYNC_DOWNLOAD as DOWNLOAD_ACCESSION2TAXID_DIR } from "../utilities/download.nf" addParams(rsync_storedir: params.accession2taxid_dir)

include { SINGLE_WGET_DOWNLOAD as DOWNLOAD_BACTERIA_SINGLE_COPY } from "../utilities/download.nf" addParams(rsync_storedir: params.bacteria_single_copy_dir)
include { SINGLE_WGET_DOWNLOAD as DOWNLOAD_BACTERIA_SINGLE_COPY_CUTOFFS } from "../utilities/download.nf" addParams(rsync_storedir: params.bacteria_single_copy_cutoffs_dir)
include { SINGLE_WGET_DOWNLOAD as DOWNLOAD_ARCHAEA_SINGLE_COPY } from "../utilities/download.nf" addParams(rsync_storedir: params.archaea_single_copy_dir)
include { SINGLE_WGET_DOWNLOAD as DOWNLOAD_ARCHAEA_SINGLE_COPY_CUTOFFS } from "../utilities/download.nf" addParams(rsync_storedir: params.archaea_single_copy_cutoffs_dir)



process DOWNLOAD_NR{
  tag "Downloading: ${filename}"
  storeDir params.nr_dir
  cpus = 1
  output:
    path "nr.dmnd"
  """
  # ascp -T -l1200m -k1 -i /home/chase/miniconda3/envs/nf-core/etc/asperaweb_id_dsa.openssh anonftp@ftp.ncbi.nlm.nih.gov:blast/db/FASTA/nr.gz /home/chase/temp/Autometa/nextflow/large_db/temp
  # docker run -v /home/chase/temp/Autometa/nextflow/large_db:/autometa/temp -it jasonkwan/autometa:dev  
  # diamond makedb --in "nr.gz" --db "nr" -p 10

  """

}

workflow DOWNLOAD_DATABASES{
    DOWNLOAD_TAXDUMP(
        "rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        "taxdump.tar.gz"
    )
    DOWNLOAD_ACCESSION2TAXID_DIR(
        "rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz",
        "prot.accession2taxid.gz"
    )
    DOWNLOAD_NR()
    
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

emit:
    ch_taxdump = DOWNLOAD_TAXDUMP.out
    ch_accession2taxid = DOWNLOAD_ACCESSION2TAXID_DIR.out
   // ch_nr = DOWNLOAD_NR.out
    ch_bacteria_single_copy = DOWNLOAD_BACTERIA_SINGLE_COPY.out
    ch_bacteria_single_copy_cutoffs = DOWNLOAD_BACTERIA_SINGLE_COPY_CUTOFFS.out
    ch_archaea_single_copy = DOWNLOAD_ARCHAEA_SINGLE_COPY.out
    ch_archaea_single_copy_cutoffs = DOWNLOAD_ARCHAEA_SINGLE_COPY_CUTOFFS.out
}

workflow {
    DOWNLOAD_DATABASES()
}