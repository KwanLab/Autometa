
process DIAMOND_MAKEDB {
    tag ' Preparing Diamond database'
    label 'process_high'

    conda "bioconda::diamond=2.0.9"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0"
    } else {
        container "quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0"
    }

    input:
        path(fasta)
        val(dbname)

    output:
        path("*.dmnd"), emit: diamond_db
        path  "versions.yml"         , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args   ?: ''
        """
        diamond makedb --in ${fasta} \\
            $args \\
            --threads ${task.cpus} \\
            --db ${dbname}


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
        END_VERSIONS
        """

}
