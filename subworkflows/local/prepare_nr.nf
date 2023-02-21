
include { DIAMOND_MAKEDB } from './../../modules/local/diamond_makedb.nf'

process DOWNLOAD_NR {
    tag "Downloading nr.gz (>100GB download. May take some time.)"
    label 'process_low'
    storeDir "${params.nr_dmnd_dir}"

    conda (params.enable_conda ? "conda-forge::rsync=3.2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path("nr.gz"), emit: singlefile

    script:
        """
        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz' 'nr.gz'

        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5' 'nr.gz.md5'

        md5sum -c *.md5
        """
}

process TEST_DOWNLOAD {
    // For development work so you don't download the entire nr.gz database
    tag "Downloading first 10,000 lines of nr.gz"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::rsync=3.2.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    output:
        path("nr.gz"), emit: singlefile

    script:
        """
        # https://github.com/nextflow-io/nextflow/issues/1564
        trap 'echo OK; exit 0;' EXIT
        curl -s ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz | zcat | head -n 10000 | gzip > nr.gz
        """
}

workflow PREPARE_NR_DB {

    main:
        if (file("${params.nr_dmnd_dir}/nr.dmnd").exists()){
            // skip huge download and db creation if nr.dmnd already exists
            out_ch = file("${params.nr_dmnd_dir}/nr.dmnd")
        } else if (file("${params.nr_dmnd_dir}/nr.gz").exists()){
            // skip huge download if nr.gz already exists
            DIAMOND_MAKEDB(file("${params.nr_dmnd_dir}/nr.gz"), "nr")
            DIAMOND_MAKEDB.out.diamond_db
                .set{out_ch}
        } else if (params.debug){
            TEST_DOWNLOAD().singlefile
                .set{nr_db_ch}
            DIAMOND_MAKEDB(nr_db_ch, "nr")
            DIAMOND_MAKEDB.out.diamond_db
                .set{out_ch}
        } else {
            DOWNLOAD_NR().singlefile
                .set{nr_db_ch}
            DIAMOND_MAKEDB(nr_db_ch, "nr")
            DIAMOND_MAKEDB.out.diamond_db
                .set{out_ch}
        }

    emit:
        diamond_db = out_ch
}
