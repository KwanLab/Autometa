

// TODO: For faster results/less I/O this could be replaced with hmmsearch
process MARKERS {
    tag "Finding markers for ${meta.id}"
    label "process_medium"

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "autometa" : null)
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
        path  '*.version.txt'               , emit: version

    script:
        def software = getSoftwareName(task.process)
        if (params.enable_conda)
        """
        exit 1
        """
        else
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

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
