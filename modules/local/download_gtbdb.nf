process DOWNLOAD_GTDB {
    tag "Downloading prot.accession2taxid.gz"
    label 'process_low'
    storeDir "${params.prot_accession2taxid_gz_dir}"

     conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path "prot.accession2taxid.gz" , emit: accession2taxid
        path "versions.yml"            , emit: versions
    when:
        task.ext.when == null || task.ext.when

    script:
        """
        autometa-config \\
            --section gtdb --option release \\
            --value ${params.gtdb_version}

        autometa-update-databases --update-gtdb

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}


workflow {
    DOWNLOAD_GTDB()
}
