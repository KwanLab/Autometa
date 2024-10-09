process RECRUIT {
    tag "sample:${meta.id}, classifier:${params.classification_method}, kmer dims:${params.classification_kmer_pca_dimensions}"
    label 'process_high'

    conda "autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    // All contigs in binning have already assigned a MAG prior to recruitment
    errorStrategy { task.exitStatus in 204 ? 'ignore' : 'terminate' }

    input:
        tuple val(meta), path(kmers), path(coverage), path(binning), path(markers), path(taxonomy)

    output:
        tuple val(meta), path("${params.kingdom}.recruitment.tsv.gz")         , emit: binning, optional: true
        tuple val(meta), path("${params.kingdom}.recruitment.main.tsv.gz")    , emit: main, optional: true
        tuple val(meta), path("${params.kingdom}.recruitment.features.tsv.gz"), emit: features, optional: true
        path  'versions.yml'                                                 , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
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
            --output-binning ${prefix}.${params.kingdom}.recruitment.tsv.gz \\
            --output-main ${prefix}.${params.kingdom}.recruitment.main.tsv.gz \\
            --output-features ${prefix}.${params.kingdom}.recruitment.features.tsv.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
        else
        """
        autometa-unclustered-recruitment \\
            --classifier ${prefix}.${params.classification_method} \\
            --kmer-dimensions ${prefix}.${params.classification_kmer_pca_dimensions} \\
            --seed 42 \\
            --taxonomy $taxonomy \\
            --kmers $kmers \\
            --coverage $coverage \\
            --binning $binning \\
            --markers $markers \\
            --output-binning ${prefix}.${params.kingdom}.recruitment.tsv.gz \\
            --output-main ${prefix}.${params.kingdom}.recruitment.main.tsv.gz \\
            --output-features ${prefix}.${params.kingdom}.recruitment.features.tsv.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
