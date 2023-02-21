
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_low'

    conda "conda-forge::python=3.8.3"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
        path samplesheet

    output:
        path '*.csv'

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        check_samplesheet.py $samplesheet samplesheet.valid.csv
        """
}
