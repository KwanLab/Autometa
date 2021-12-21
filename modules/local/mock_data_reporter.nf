// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MOCK_DATA_REPORT {

    tag ' Preparing mock data report'
    label 'process_low'

    publishDir "${params.outdir_internal}/${options.publish_dir}"

    // container  TODO: The R "Rocker" Docker images don't have ps which is required by Nextflow so this may have to be a custom image?
    // docker build https://gist.githubusercontent.com/chasemc/818111640daae05beb2b070641aa33fb/raw/09107704fa60a6311fb09542c8b99b848e168ea3/Dockerfile --tag mock_data_reporter


    input:
        tuple val(meta), path(bins_path), path(assembly_to_locus_path), path(assembly_report_path)
        path(rmarkdown_file)

    output:
        tuple val(meta), path("*.html"), emit: results


    script:
        """
        #!/usr/bin/env Rscript

        packages <- c("markdown","data.table", "ggplot2", "plotly", "crosstalk", "magrittr", "DT", "stringi")

        for (i in packages) {
          if (!requireNamespace(i)) {
            install.packages(i)
          }
          library(i, character.only = T)
        }

        rmarkdown::render("${rmarkdown_file}",
                        params = list(
                          bins_path= "${bins_path}",
                          assembly_to_locus_path = "${assembly_to_locus_path}",
                          assembly_report_path = "${assembly_report_path}",
                          genus=FALSE
                        ),
                        knit_root_dir=getwd(),
                        output_dir=getwd(),
                        output_file="mock_data_report_by_assembly.html")

        rmarkdown::render("${rmarkdown_file}",
                        params = list(
                          bins_path= "${bins_path}",
                          assembly_to_locus_path = "${assembly_to_locus_path}",
                          assembly_report_path = "${assembly_report_path}",
                          genus=TRUE
                        ),
                        knit_root_dir=getwd(),
                        output_dir=getwd(),
                        output_file="mock_data_report_by_genus.html")

        """
}
