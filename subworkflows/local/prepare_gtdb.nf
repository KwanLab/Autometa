
include { DIAMOND_MAKEDB as GTDB_MAKEDB } from './../../modules/local/diamond_makedb.nf'


process DOWNLOAD_GTDB {
    tag "Downloading GTDB database version ${params.gtdb_version}"
    label 'process_low'
    storeDir "${params.gtdb_dir}"
    
    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path 'autometa_formatted_gtdb-version-*.faa.gz' , emit: gtdb_formated_faa
        path 'gtdb_taxdump-version-*/*'                 , emit: gtdb_taxdump_directory
        path "versions.yml"                             , emit: versions

    script:
        """
        autometa-download-gtdb --version $params.gtdb_version --outdir '.'

        rm gtdb-taxdump-version-*.tar.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
            gtdb: $params.gtdb_version
        END_VERSIONS
        """
}

workflow PREPARE_GTDB_DB {

    main:
        ch_versions = Channel.empty()

        // if use_gtdb and large_downloads_permission is set to true, download the GTDB database
        if (params.use_gtdb && params.large_downloads_permission) {

            DOWNLOAD_GTDB()
            ch_versions = ch_versions.mix(DOWNLOAD_GTDB.out.versions)

            // get the single gtdb_formated_faa file and create the string e.g. autometa_formatted_gtdb-version-220.db from  autometa_formatted_gtdb-version-220.faa.gz
            dbname = DOWNLOAD_GTDB.out.gtdb_formated_faa.getName().replaceFirst(/\.gz$/, '').replaceFirst(/\.faa$/, '.dmnd')           
            GTDB_MAKEDB(DOWNLOAD_GTDB.out.gtdb_formated_faa, dbname)
            ch_versions = ch_versions.mix(GTDB_MAKEDB.out.versions)
            
        } else {
            println '\033[0;34m `--large_downloads_permission` is set to false. Skipping GTDB database download. \033[0m'
        }

    emit:
        diamond_db = GTDB_MAKEDB.out.diamond_db
        gtdb_taxdump_directory = DOWNLOAD_GTDB.out.gtdb_taxdump_directory
        versions = ch_versions
}
