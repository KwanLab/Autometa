
process COUNT_KMERS {
    tag "Counting ${params.kmer_size}-mers for ${meta.id}"
    label 'process_medium'
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("kmers.tsv")           , emit: counts
        path  '*.version.txt'                        , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-kmers \\
            --fasta $metagenome \\
            --kmers "kmers.tsv" \\
            --size "${params.kmer_size}" \\
            --cpus "${task.cpus}" \\
            --seed 42

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
