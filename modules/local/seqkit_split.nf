process SEQKIT_SPLIT {
    tag "Splitting $meta.id for parallel processing"
    label 'process_medium'

    conda "bioconda::seqkit=0.16.1"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0"
    }

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path("outfolder/*")    , emit: fasta
        path "*.version.txt"                    , emit: version

    when:
        task.ext.when == null || task.ext.when

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        """
        seqkit \\
            split \\
            ${fasta} \\
            ${options.args} \\
            ${options.args2} \\
            -O outfolder

        seqkit version | sed 's/seqkit v//g' > ${software}.version.txt
        """
}
