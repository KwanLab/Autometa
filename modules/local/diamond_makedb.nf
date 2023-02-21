
process DIAMOND_MAKEDB {
    tag ' Preparing Diamond database'
    label 'process_high'

    storeDir "${params.nr_dmnd_dir}"

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
        path  "*.version.txt"         , emit: version

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        diamond makedb --in ${fasta} \\
            $options.args \\
            --threads ${task.cpus} \\
            --db ${dbname}

        diamond version | sed 's/^.*diamond version //' > diamond.version.txt
        """
}
