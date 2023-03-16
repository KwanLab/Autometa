// TODO: For faster results/less I/O this could be replaced with hmmsearch
process MARKERS {
    tag "Finding markers for ${meta.id}, taxon:${meta.taxon}"
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
        tuple val(meta), path("*.markers.tsv.gz")  , emit: markers_tsv
        tuple val(meta), path("*.hmmscan.tsv.gz")  , emit: hmmscan_tsv
        path  'versions.yml'                    , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}.${meta.taxon}"
        """
        autometa-markers \\
            --orfs $orfs \\
            --hmmscan ${prefix}.hmmscan.tsv \\
            --out ${prefix}.markers.tsv \\
            --kingdom ${meta.taxon} \\
            --parallel \\
            --cpus ${task.cpus} \\
            --seed 42 \\
            --hmmdb "/scratch/dbs/markers/${meta.taxon}.single_copy.hmm" \\
            --cutoffs "/scratch/dbs/markers/${meta.taxon}.single_copy.cutoffs"


        gzip -6  ${prefix}.hmmscan.tsv
        gzip -6  ${prefix}.markers.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
