process DOWNLOAD_ACESSION2TAXID {
    tag "Downloading prot.accession2taxid.gz"
    label 'process_low'
    storeDir "${params.prot_accession2taxid_gz_dir}"

    conda "conda-forge::rsync=3.2.3"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        // hack nf-core options.args3 and use for output name
        path "prot.accession2taxid.gz" , emit: accession2taxid
        path "versions.yml"            , emit: versions
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz' 'prot.accession2taxid.gz'

        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz.md5' 'prot.accession2taxid.gz.md5'

        md5sum -c *.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rsync: \$(rsync --version | head -n1 | sed 's/^rsync  version //' | sed 's/\s.*//')
        END_VERSIONS
        """
    stub:
        """
        # https://github.com/nextflow-io/nextflow/issues/1564
        trap 'echo OK; exit 0;' EXIT
        curl -s ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz | zcat | head -n 1000 |gzip > prot.accession2taxid.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rsync: \$(rsync --version | head -n1 | sed 's/^rsync  version //' | sed 's/\s.*//')
        END_VERSIONS
        """
}

process DOWNLOAD_TAXDUMP {
    tag "Downloading taxdump.tar.gz"
    label 'process_low'

    conda "conda-forge::rsync=3.2.3"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path "*.dmp"        , emit: taxdump_files
        path "versions.yml" , emit: versions

    when:
        task.ext.when == null || task.ext.when

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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rsync: \$(rsync --version | head -n1 | sed 's/^rsync  version //' | sed 's/\s.*//')
        END_VERSIONS
        """
}

workflow PREPARE_TAXONOMY_DATABASES {
    main:
        ch_versions = Channel.empty()

        taxdump_dir = file(params.taxdump_tar_gz_dir)
        taxdump_dir_files = taxdump_dir.list()
        expected_files = ['citations.dmp', 'delnodes.dmp', 'division.dmp', 'gencode.dmp', 'merged.dmp', 'names.dmp', 'nodes.dmp']

        dmp_files = file("${params.taxdump_tar_gz_dir}/*.dmp")
        taxonomy_files_exist = dmp_files.name.containsAll(expected_files)

        if (taxonomy_files_exist){
            taxdump_files = dmp_files
        } else {
            DOWNLOAD_TAXDUMP()
            ch_versions = ch_versions.mix(DOWNLOAD_TAXDUMP.out.versions)
            DOWNLOAD_TAXDUMP.out.taxdump_files
                .set{taxdump_files}
        }

        expected_files2 = ['prot.accession2taxid.gz']

        taxonomy_files_exist2 = file("${params.prot_accession2taxid_gz_dir}/*.dmp").name.containsAll(expected_files2)

        if (taxonomy_files_exist2){
            prot_accession2taxid_ch = accession2taxid_dir_files
        } else {
            DOWNLOAD_ACESSION2TAXID()
            DOWNLOAD_ACESSION2TAXID.out.accession2taxid
                .set{prot_accession2taxid_ch}
            ch_versions = ch_versions.mix(DOWNLOAD_ACESSION2TAXID.out.versions)

        }

    emit:
        taxdump_files = taxdump_files
        prot_accession2taxid = prot_accession2taxid_ch
        versions = ch_versions

}

