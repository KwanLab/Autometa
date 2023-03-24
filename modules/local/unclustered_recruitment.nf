process UNCLUSTERED_RECRUIT {
    tag "sample:${meta.id}, classifier:${params.classification_method}, kmer dims:${params.classification_kmer_pca_dimensions}"
    label 'process_high'

    conda "bioconda::autometa>=${params.autometa_image_tag}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    // All contigs in binning have already assigned a MAG prior to recruitment
    errorStrategy { task.exitStatus in 204 ? 'ignore' : 'terminate' }

    input:
        tuple val(meta), path(kmers), path(coverage), path(binning), path(markers), path(taxonomy)

    output:
        tuple val(meta), path("*.recruitment.tsv.gz")           , emit: binning   , optional: true
        tuple val(meta), path("*.recruitment.main.tsv.gz")      , emit: main      , optional: true
        tuple val(meta), path("*.recruitment.features.tsv.gz")  , emit: features  , optional: true
        path  'versions.yml'                                    , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}.${meta.taxon}"
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
                --output-binning ${prefix}.recruitment.tsv.gz \\
                --output-main ${prefix}.recruitment.main.tsv.gz \\
                --output-features ${prefix}.recruitment.features.tsv.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                autometa: \$(autometa --version | sed -e 's/autometa: //g')
            END_VERSIONS
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
                --output-binning ${prefix}.recruitment.tsv.gz \\
                --output-main ${prefix}.recruitment.main.tsv.gz \\
                --output-features ${prefix}.recruitment.features.tsv.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                autometa: \$(autometa --version | sed -e 's/autometa: //g')
            END_VERSIONS
            """
}
