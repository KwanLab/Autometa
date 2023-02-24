process MOCK_DATA_REPORT {

    tag 'Preparing mock data report'
    label 'process_low'

    publishDir "${options.publish_dir}", mode: params.publish_dir_mode

    container "jasonkwan/autometa-nf-modules-mock_data_reporter:main"

    input:
        tuple val(meta), path(bins_path), path(assembly_to_locus_path), path(assembly_report_path)
        path(rmarkdown_file)

    output:
        tuple val(meta), path("*.html") , emit: results
        path "versions.yml"             , emit: versions


    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mock_data_report.R ${rmarkdown_file} ${bins_path} ${assembly_to_locus_path} ${assembly_report_path}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            R: 'For R and packages, see docker: jasonkwan/autometa-nf-modules-mock_data_reporter:main'
        END_VERSIONS

        """
}
