
/*
=======================
TODO: Not yet implemented
"Cutoffs" need to be downloaded/provided to this process
=======================
*/


process HMMER_HMMSEARCH_FILTER {
    tag "Filtering marker hmms in ${meta.id}"
    label 'process_medium'

    // if ( params.num_splits < 2 ) {
    // if running in parallel, the results are published from the process
    // that merges the individual results from this process

    conda "autometa"
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
        path "*.version.txt"                           , emit: version

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

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
