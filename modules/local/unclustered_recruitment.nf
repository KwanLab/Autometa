// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process RECRUIT {
    tag "sample:${meta.id}, classifier:${params.classification_method}, kmer dims:${params.classification_kmer_pca_dimensions}"
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
        tuple val(meta), path ("${params.kingdom}.recruitment.tsv.gz")         , emit: binning, optional: true
        tuple val(meta), path ("${params.kingdom}.recruitment.main.tsv.gz")    , emit: main, optional: true
        tuple val(meta), path ("${params.kingdom}.recruitment.features.tsv.gz"), emit: features, optional: true
        path  '*.version.txt'                                                  , emit: version

    script:
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
            --output-main ${params.kingdom}.recruitment.main.tsv.gz \\
            --output-features ${params.kingdom}.recruitment.features.tsv.gz

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
            --output-main ${params.kingdom}.recruitment.main.tsv.gz \\
            --output-features ${params.kingdom}.recruitment.features.tsv.gz

        echo "TODO" > autometa.version.txt
        """
}
