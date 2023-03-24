
/*
=======================
TODO: Not yet implemented
This should speed up hmm searches, however as of now 2 things are needed:
1: cutoff values need to be downloaded/provided to the next process that reads
the results of this process
2: The cutoffs would need to be determined again using the -Z flag of hmmsearch
=======================


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
        path "versions.yml"                , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmsearch -h | grep -o '^# HMMER [0-9.]*' | sed 's/^# HMMER *//')
    END_VERSIONS
    """
}


process HMMER_HMMSEARCH_FILTER {
    tag "Filtering marker hmms in ${meta.id}"
    label 'process_medium'

    // if ( params.num_splits < 2 ) {
    // if running in parallel, the results are published from the process
    // that merges the individual results from this process

    conda "bioconda::autometa=${params.autometa_image_tag}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(domtblout), path(fasta)
        path("bacteria.single_copy.cutoffs") // FIXME: This should take cutoffs corresponding to domtblout

    output:
        tuple val(meta), path("markers.tsv"), emit: markers_tsv
        path "versions.yml"                           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        autometa-hmmsearch-filter \\
            --domtblout "$domtblout" \\
            --cutoffs  TODO:"Cutoffs" need to be downloaded/provided to this process \\
            --seqdb "$fasta" \\
            --out "markers.tsv"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
*/
