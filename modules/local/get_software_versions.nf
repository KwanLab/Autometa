// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process GET_SOFTWARE_VERSIONS {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    cache false

    input:
        path versions

    output:
        path "software_versions.tsv"     , emit: tsv
        path 'software_versions_mqc.yaml', emit: yaml
        path  '*.version.txt'            , emit: version

    script:
        // Add soft-links to original FastQs for consistent naming in pipeline
        def software = getSoftwareName(task.process)
        """
        echo $workflow.manifest.version > pipeline.version.txt
        echo $workflow.nextflow.version > nextflow.version.txt
        scrape_software_versions.py &> software_versions_mqc.yaml

        echo "make linter happy" > autometa.version.txt
        """
}
