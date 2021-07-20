// this file probably needs to be reevaluated, but from a python-first
// perspective since the python code assumes file/directory structure
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.taxdump_tar_gz_dir = [:]
params.prot_accession2taxid_gz_dir = [:]
options = initOptions(params.options)




process TEST_DOWNLOAD {
    // For development work so you don't download the entire nr.gz database
    tag "Downloading first 10,000 lines of nr.gz"
    label 'process_low'
    storeDir "${params.prot_accession2taxid_gz_dir}"

    conda (params.enable_conda ? "conda-forge::rsync=3.2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
         container "jason-c-kwan/autometa:nfcore"
    }

    output:
    path("prot.accession2taxid"), emit: singlefile

    """
    # https://github.com/nextflow-io/nextflow/issues/1564
    trap 'echo OK; exit 0;' EXIT
    curl -s ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz | zcat | head -n 1000 > prot.accession2taxid
    """
}

process DOWNLOAD_ACESSION2TAXID {
    tag "Downloading prot.accession2taxid.gz"
    label 'process_low'
    storeDir "${params.prot_accession2taxid_gz_dir}"

    conda (params.enable_conda ? "conda-forge::rsync=3.2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
         container "jason-c-kwan/autometa:nfcore"
    }

    output:
    // hack nf-core options.args3 and use for output name
    path "prot.accession2taxid" , emit: singlefile
    path  "*.version.txt"   , emit: version

    """
    rsync -a \
        --quiet \
        'rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz' 'prot.accession2taxid.gz'
        'rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5' 'prot.accession2taxid.gz.md5'

    md5sum -c *.md5
    gunzip prot.accession2taxid.gz
    rsync --version | head -n1 > rsync.version.txt
    """
}


process DOWNLOAD_TAXDUMP {
    tag "Downloading taxdump.tar.gz"
    label 'process_low'
    storeDir "${params.taxdump_tar_gz_dir}"

    conda (params.enable_conda ? "conda-forge::rsync=3.2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
         container "jason-c-kwan/autometa:nfcore"
    }

    output:
    // hack nf-core options.args3 and use for output name
    path "*" , emit: singlefile
    path  "*.version.txt"   , emit: version

    """
    rsync -a \
        --quiet \
        'rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz' 'taxdump.tar.gz'

    rsync -a \
        --quiet \
        'rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5' 'taxdump.tar.gz.md5'

    md5sum -c *.md5
   rm 'taxdump.tar.gz.md5'
   tar -xf taxdump.tar.gz
   rm taxdump.tar.gz

    rsync --version | head -n1 > rsync.version.txt
    """
}


workflow PREPARE_TAXONOMY_DATABASES {
    main:

         if (params.debug){
            TEST_DOWNLOAD().singlefile
                .set{prot_accession2taxid_ch}
        }
        else {
            DOWNLOAD_ACESSION2TAXID().singlefile
                .set{prot_accession2taxid_ch}
        }
        DOWNLOAD_TAXDUMP()

    emit:
        taxdump = DOWNLOAD_TAXDUMP.out.singlefile
        prot_accession2taxid = prot_accession2taxid_ch

}

