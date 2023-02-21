
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

    conda "bioconda::hmmer=3.3.2"
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

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        # hmmsearch can'ts use or pipe in gzipped fasta

        zcat "${fasta}" > temp.fa

        hmmsearch \\
            --domtblout "${hmm.simpleName}.domtblout" \\
            --cpu $task.cpus \\
            $args \\
            $args2 \\
            "${hmm}" \\
            temp.fa > /dev/null 2>&1

        echo \$(hmmalign -h | grep -o '^# HMMER [0-9.]*') | sed 's/^# HMMER *//' > HMMER.version.txt
        """
}
