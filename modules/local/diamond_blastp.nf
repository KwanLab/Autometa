process DIAMOND_BLASTP {
    tag "Aligning ORFS in ${meta.id} against ${diamond_database}"
    label 'process_high'
    // Old diamond manual suggested *NOT* running in parallel... so we are setting maxForks to 1 here.
    // TODO: There appears to be features for multiprocessing available now
    // See: https://github.com/bbuchfink/diamond/wiki/6.-Distributed-computing
    maxForks 1

    conda (params.enable_conda ? "bioconda::diamond=2.0.14" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/diamond:2.0.14--hdcc8f71_0"
    } else {
        container "quay.io/biocontainers/diamond:2.0.14--hdcc8f71_0"
    }

    input:
        tuple val(meta), path(protein_fasta)
        path(diamond_database)

    output:
        tuple val(meta), path("blastp.tsv"), emit: diamond_results
        path "*.version.txt"               , emit: version

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        """
        diamond blastp $args \\
            --query ${protein_fasta} \\
            --db ${diamond_database} \\
            --threads ${task.cpus} \\
            --out blastp.tsv

        diamond version | sed 's/^.*diamond version //' > diamond.version.txt
        """
}
