// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BIN_CONTIGS {
    tag "Performing Autometa binning on ${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(kmers), path(coverage), path(gc_content), path(markers), path(taxonomy)

    output:
        tuple val(meta), path("${params.kingdom}.binning.tsv.gz"), emit: binning
        tuple val(meta), path("${params.kingdom}.main.tsv.gz")   , emit: main
        path  '*.version.txt'                                               , emit: version

    script:
        def software = getSoftwareName(task.process)
        taxonomy_call = params.taxonomy_aware ? "--taxonomy $taxonomy" : "" // https://github.com/nextflow-io/nextflow/issues/1694#issuecomment-683272275
        """
        autometa-binning \\
            --kmers $kmers \\
            --coverages $coverage \\
            --gc-content $gc_content \\
            --markers $markers \\
            --output-binning ${params.kingdom}.binning.tsv.gz \\
            --output-main ${params.kingdom}.main.tsv.gz \\
            --clustering-method ${params.clustering_method} \\
            --completeness ${params.completeness} \\
            --purity ${params.purity} \\
            $taxonomy_call \\
            --cov-stddev-limit ${params.cov_stddev_limit} \\
            --gc-stddev-limit ${params.gc_stddev_limit} \\
            --starting-rank ${params.binning_starting_rank} \\
            --domain ${params.kingdom}

        echo "TODO" > autometa.version.txt
        """
}
