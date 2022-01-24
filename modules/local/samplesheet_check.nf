
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'
    publishDir "${params.outdir}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
        path samplesheet

    output:
        path '*.csv'

    script:
        """
        check_samplesheet.py $samplesheet samplesheet.valid.csv
        """
}
