
/*
=======================
TODO: Not yet implemented
This should speed up hmm searches, however as of now 2 things are needed:
1: cutoff values need to be downloaded/provided to the next process that reads
the results of this process
2: The cutoffs would need to be determined again using the -Z flag of hmmsearch
=======================
*/

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
        path "*.version.txt"                , emit: version

    script:
        def software = getSoftwareName(task.process)
        def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
        def fastacmd = fasta.getExtension() == 'gz' ? "gunzip -c $fasta" : "cat $fasta"
        """
        hmmsearch \\
            --domtblout "${hmm.simpleName}.domtblout" \\
            ${options.args} \\
            ${options.args2} \\
            $hmm \\
            $fasta > /dev/null 2>&1

        echo \$(hmmalign -h | grep -o '^# HMMER [0-9.]*') | sed 's/^# HMMER *//' > HMMER.version.txt
        """
}
