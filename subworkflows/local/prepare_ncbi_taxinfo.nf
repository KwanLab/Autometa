// this file probably needs to be reevaluated, but from a python-first
// perspective since the python code assumes file/directory structure
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
params.taxdump_tar_gz_dir = [:]
params.prot_accession2taxid_gz_dir = [:]
options = initOptions(params.options)

process TEST_DOWNLOAD {
    // For development work so you don't download the entire prot.accession2taxid.gz database
    tag "Downloading first 10,000 lines of prot.accession2taxid.gz"
    label 'process_low'
    storeDir "${params.prot_accession2taxid_gz_dir}"

    conda (params.enable_conda ? "conda-forge::rsync=3.2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path("prot.accession2taxid"), emit: singlefile

    script:
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
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        // hack nf-core options.args3 and use for output name
        path "prot.accession2taxid.gz" , emit: accession2taxid
        path  "*.version.txt"   , emit: version
    script:
        """
        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz' 'prot.accession2taxid.gz'

        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5' 'prot.accession2taxid.gz.md5'

        md5sum -c *.md5

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
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path "*" , emit: taxdump_files
        path  "*.version.txt"   , emit: version

    script:
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
        expected_files = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'names.dmp', 'nodes.dmp']
        taxdump_dir = file(params.taxdump_tar_gz_dir)
        taxdump_dir_files = []
        taxdump_dir.eachFile { item ->
            if( item.isFile() ) {
              //  println "${item.getName()}"
                taxdump_dir_files.add(item.getName())
            }
            else if( item.isDirectory() ) {
                //println "${item.getName()} - DIR"
            }
        }
        if (taxdump_dir_files.containsAll(expected_files)){
            taxdump_files = taxdump_dir_files
        } else {
            DOWNLOAD_TAXDUMP()
            DOWNLOAD_TAXDUMP.out.taxdump_files
                .set{taxdump_files}
        }

        accession2taxid_dir = file(params.prot_accession2taxid_gz_dir)
        accession2taxid_dir_files = []
        expected_files = ['prot.accession2taxid']

        accession2taxid_dir.eachFile { item ->
            if( item.isFile() ) {
              //  println "${item.getName()}"
                accession2taxid_dir_files.add(item.getName())
            }
            else if( item.isDirectory() ) {
                //println "${item.getName()} - DIR"
            }
        }
        if (accession2taxid_dir_files.containsAll(expected_files)){
            prot_accession2taxid_ch = accession2taxid_dir_files
        } else if (params.debug){
            TEST_DOWNLOAD().singlefile
                .set{prot_accession2taxid_ch}
        } else {
            DOWNLOAD_ACESSION2TAXID().accession2taxid
                .set{prot_accession2taxid_ch}
        }

    emit:
        taxdump = taxdump_files
        prot_accession2taxid = prot_accession2taxid_ch

}

