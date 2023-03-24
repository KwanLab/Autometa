process BINNING {
    tag "sample:${meta.id}, taxon:${meta.taxon}, clustering:${params.clustering_method}, completeness:${params.completeness}, purity:${params.purity}, cov.std.dev.:${params.cov_stddev_limit}, gc.std.dev.:${params.gc_stddev_limit}"
    label 'process_high'

    conda "bioconda::autometa>=${params.autometa_image_tag}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    // No markers were annotated for contigs in the table
    errorStrategy { task.exitStatus in 204 ? 'ignore' : 'terminate' }

    input:
        tuple val(meta), path(kmers), path(coverage), path(gc_content), path(markers), path(taxonomy)

    output:
        tuple val(meta), path("*.binning.tsv.gz")        , emit: binning
        tuple val(meta), path("*.binning.main.tsv.gz")   , emit: main
        path 'versions.yml'                              , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        taxonomy_call = params.taxonomy_aware ? "--taxonomy $taxonomy" : "" // https://github.com/nextflow-io/nextflow/issues/1694#issuecomment-683272275
        def prefix = task.ext.prefix ?: "${meta.id}.${meta.taxon}"
        """
        autometa-binning \\
            --kmers $kmers \\
            --coverages $coverage \\
            --gc-content $gc_content \\
            --markers ${markers} \\
            --output-binning ${prefix}.binning.tsv.gz \\
            --output-main ${prefix}.binning.main.tsv.gz \\
            --clustering-method ${params.clustering_method} \\
            --completeness ${params.completeness} \\
            --purity ${params.purity} \\
            $taxonomy_call \\
            --cov-stddev-limit ${params.cov_stddev_limit} \\
            --gc-stddev-limit ${params.gc_stddev_limit} \\
            --starting-rank ${params.binning_starting_rank} \\
            --cpus ${task.cpus} \\
            --rank-filter superkingdom \\
            --rank-name-filter ${meta.taxon} \\
            --verbose

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
