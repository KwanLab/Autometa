// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process HMMER_HMMSEARCH {
    tag "Annotating ORFs in $meta.id"
    label 'process_medium'

    // no publishdir

    conda (params.enable_conda ? "bioconda::hmmer=3.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/hmmer:3.3.2--h1b792b2_1"
    } else {
        container "quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1"
    }

    input:
        tuple val(meta), path(fasta)
        path(hmm)

    output:
        tuple val(meta), path("*.domtblout"), emit: domtblout
        path "*.version.txt"               , emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        def fastacmd = fasta.getExtension() == 'gz' ? "gunzip -c $fasta" : "cat $fasta"
        """
        hmmsearch \\
            --domtblout "${meta.id}.domtblout" \\
            --tblout "${meta.id}.tblout" \\
            --pfamtblout "${meta.id}.pfamtblout" \\
            ${options.args} \\
            ${options.args2} \\
            $hmm \\
            $fasta > /dev/null 2>&1

        echo \$(hmmalign -h | grep -o '^# HMMER [0-9.]*') | sed 's/^# HMMER *//' > HMMER.version.txt
        """
}
