// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RECRUIT {
    tag "Performing Autometa unclustered recruitment on ${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(kmers), path(coverage), path(binning), path(markers), path(taxonomy)

    output:
        tuple val(meta), path ("${params.kingdom}.recruitment.tsv.gz")     , emit: binning, optional: true
        tuple val(meta), path ("${params.kingdom}.recruitment.main.tsv.gz"), emit: main, optional: true
        path  '*.version.txt'                                              , emit: version

    script:
        // Add soft-links to original FastQs for consistent naming in pipeline
        def software = getSoftwareName(task.process)
        if (!params.taxonomy_aware)
        """
        autometa-unclustered-recruitment \\
            --classifier ${params.classification_method} \\
            --kmer-dimensions ${params.classification_kmer_pca_dimensions} \\
            --seed 42 \\
            --kmers $kmers \\
            --coverage $coverage \\
            --binning $binning \\
            --markers $markers \\
            --output-binning ${params.kingdom}.recruitment.tsv.gz \\
            --output-main ${params.kingdom}.recruitment.main.tsv.gz
        echo "TODO" > autometa.version.txt
        """
        else
        """
        autometa-unclustered-recruitment \\
            --classifier ${params.classification_method} \\
            --kmer-dimensions ${params.classification_kmer_pca_dimensions} \\
            --seed 42 \\
            --taxonomy $taxonomy \\
            --kmers $kmers \\
            --coverage $coverage \\
            --binning $binning \\
            --markers $markers \\
            --output-binning ${params.kingdom}.recruitment.tsv.gz \\
            --output-main ${params.kingdom}.recruitment.main.tsv.gz

        echo "TODO" > autometa.version.txt
        """
}
