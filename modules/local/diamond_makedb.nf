// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'


params.nr_dmnd_dir = nullprocess DIAMOND_MAKEDB {
    tag ' Preparing Diamond database'
    label 'process_high'

    storeDir "${params.nr_dmnd_dir}"

    conda (params.enable_conda ? "bioconda::diamond=2.0.9" : null)
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

    script:
        def software = getSoftwareName(task.process)
        """
        diamond makedb --in ${fasta} \\
            $options.args \\
            --threads ${task.cpus} \\
            --db ${dbname}

        diamond version | sed 's/^.*diamond version //' > diamond.version.txt
        """
}
