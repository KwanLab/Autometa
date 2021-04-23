#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: "${params.publish_dir_mode}"

    input:
    path output_docs 

    output:
    path 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}