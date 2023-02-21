process MOCK_DATA_REPORT {

    tag 'Preparing mock data report'
    label 'process_low'

    publishDir "${options.publish_dir}", mode: params.publish_dir_mode

    container "jasonkwan/autometa-nf-modules-mock_data_reporter:main"


    input:
        tuple val(meta), path(bins_path), path(assembly_to_locus_path), path(assembly_report_path)
        path(rmarkdown_file)

    output:
        tuple val(meta), path("*.html"), emit: results


    when:
        task.ext.when == null || task.ext.when

    script:
        """
        #!/usr/bin/env Rscript
        rmarkdown::render(
          input="${rmarkdown_file}",
          params=list(
            bins_path="${bins_path}",
            assembly_to_locus_path="${assembly_to_locus_path}",
            assembly_report_path="${assembly_report_path}",
            genus=FALSE
          ),
          knit_root_dir=getwd(),
          output_dir=getwd(),
          output_file="mock_data_report_by_assembly.html"
        )

        rmarkdown::render(
          input="${rmarkdown_file}",
          params=list(
            bins_path= "${bins_path}",
            assembly_to_locus_path = "${assembly_to_locus_path}",
            assembly_report_path = "${assembly_report_path}",
            genus=TRUE
          ),
          knit_root_dir=getwd(),
          output_dir=getwd(),
          output_file="mock_data_report_by_genus.html"
        )

        """
}
