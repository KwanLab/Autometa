

// TODO: For faster results/less I/O this could be replaced with hmmsearch
process MARKERS {
    tag "Finding markers for ${meta.id}"
    label "process_medium"


    conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
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
        tuple val(meta), path("${params.kingdom}.markers.tsv"), emit: markers_tsv
        tuple val(meta), path("${params.kingdom}.hmmscan.tsv"), emit: hmmscan_tsv
        path  'versions.yml'               , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:

     // TODO: nf-core linter now checks and flags any instances of 'params.enable_conda'. Is this check below still necessary?
     //   if (params.enable_conda)
     //   """
     //   exit 1
     //   """
     //   else
        """
        autometa-markers \\
            --orfs $orfs \\
            --hmmscan ${params.kingdom}.hmmscan.tsv \\
            --out ${params.kingdom}.markers.tsv \\
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
