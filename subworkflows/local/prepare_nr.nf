
include { DIAMOND_MAKEDB } from './../../modules/local/diamond_makedb.nf'

process DOWNLOAD_NR {
    tag "Downloading nr.gz (>100GB download. May take some time.)"
    label 'process_low'
    label 'process_long'

    output:
        path("nr.gz")       , emit: ncbi_nr_fasta
        path("nr.gz.md5")   , emit: md5
        path "versions.yml" , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz' 'nr.gz'

        rsync -a \\
            --quiet \\
            'rsync://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5' 'nr.gz.md5'

        md5sum -c *.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            rsync: \$(rsync --version | head -n1 | sed 's/^rsync  version //' | sed 's/\s.*//')
        END_VERSIONS
        """
    stub:
        """
        curl https://raw.githubusercontent.com/chasemc/autometa_test_data/main/example_1/nr.fa.gz  > nr.gz
        echo 'faac6c3d544e3d89f5b99377814f0c37 nr.gz' > 'nr.gz.md5'
        md5sum -c 'nr.gz.md5'
        touch versions.yml
        """
}

workflow PREPARE_NR_DB {

    main:
        ch_versions = Channel.empty()

        // TODO: this if/else can be simplified


        if (file("${params.nr_dmnd_dir}/nr.dmnd").exists()){

            // skip huge download and db creation if nr.dmnd already exists
            out_ch = file("${params.nr_dmnd_dir}/nr.dmnd")

        } else if (file("${params.nr_dmnd_dir}/nr.gz").exists()){

            // skip huge download if nr.gz already exists
            DIAMOND_MAKEDB(file("${params.nr_dmnd_dir}/nr.gz"), "nr")
            ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)

            out_ch = DIAMOND_MAKEDB.out.diamond_db

        } else if (params.large_downloads_permission || workflow.stubRun) {

            DOWNLOAD_NR()
            ch_versions = ch_versions.mix(DOWNLOAD_NR.out.versions)

            DIAMOND_MAKEDB(DOWNLOAD_NR.out.ncbi_nr_fasta, "nr")
            ch_versions = ch_versions.mix(DIAMOND_MAKEDB.out.versions)

            out_ch = DIAMOND_MAKEDB.out.diamond_db

        } else {
            println '\033[0;34m Neither nr.dmnd or nr.gz were found and `--large_downloads_permission` is set to false. \033[0m'
        }

    emit:
        diamond_db = out_ch
        versions = ch_versions
}
