

// TODO: For faster results/less I/O this could be replaced with hmmsearch
process MARKERS {
    tag "Finding markers for ${meta.id}"
    label "process_medium"


    conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    // No markers in kingdom.hmmscan.tsv pass cutoff thresholds
    errorStrategy { task.exitStatus in 204 ? 'ignore' : 'terminate' }

    input:
        tuple val(meta), path(orfs)
        //path(hmmdb) currently only inside docker
        //path(cutoffs) currently only inside docker

    output:
        tuple val(meta), path("*.markers.tsv")  , emit: markers_tsv
        tuple val(meta), path("*.hmmscan.tsv")  , emit: hmmscan_tsv
        path  'versions.yml'                    , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        autometa-markers \\
            --orfs $orfs \\
            --hmmscan ${prefix}.${params.kingdom}.hmmscan.tsv \\
            --out ${prefix}.${params.kingdom}.markers.tsv \\
            --kingdom ${params.kingdom} \\
            --parallel \\
            --cpus ${task.cpus} \\
            --seed 42 \\
            --hmmdb "/scratch/dbs/markers/${params.kingdom}.single_copy.hmm" \\
            --cutoffs "/scratch/dbs/markers/${params.kingdom}.single_copy.cutoffs"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
