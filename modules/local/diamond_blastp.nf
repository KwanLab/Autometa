process DIAMOND_BLASTP {
    tag "Aligning ORFS in ${meta.id} against ${diamond_database}"
    label 'process_high'
    // Old diamond manual suggested *NOT* running in parallel... so we are setting maxForks to 1 here.
    // TODO: There appears to be features for multiprocessing available now
    // See: https://github.com/bbuchfink/diamond/wiki/6.-Distributed-computing
    maxForks 1

    conda "bioconda::diamond=2.0.14"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/diamond:2.0.14--hdcc8f71_0"
    } else {
        container "quay.io/biocontainers/diamond:2.0.14--hdcc8f71_0"
    }

    input:
        tuple val(meta), path(protein_fasta)
        path(diamond_database)

    output:
        tuple val(meta), path("*blastp.tsv"), emit: diamond_results
        path "versions.yml"                 , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        diamond blastp $args \\
            --query ${protein_fasta} \\
            --db ${diamond_database} \\
            --threads ${task.cpus} \\
            --out ${prefix}.blastp.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
        END_VERSIONS
        """
}
