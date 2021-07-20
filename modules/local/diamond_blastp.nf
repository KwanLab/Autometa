// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DIAMOND_BLASTP {
    tag "Annotating ORFS in ${meta.id}"
    label 'process_high'
    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::diamond=2.0.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0"
    } else {
        container "quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0"
    }

    input:
        tuple val(meta), path(protein_fasta)
        path(diamond_database)

    output:
        tuple val(meta), path("${meta.id}.blastp.tsv"), emit: diamond_results
        path "*.version.txt"               , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        diamond blastp $options.args \\
            --query ${protein_fasta} \\
            --db ${diamond_database} \\
            --threads ${task.cpus} \\
            --out ${meta.id}.blastp.tsv

        diamond version | sed 's/^.*diamond version //' > diamond.version.txt
        """
}
