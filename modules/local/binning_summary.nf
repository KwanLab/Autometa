process BINNING_SUMMARY {
    tag "Gathering binning summary for ${meta.id}"
    label 'process_high'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(binning_main), path(markers), path(metagenome), val(binning_column)
        path("*")

    output:
        tuple val(meta), path("*metabin_stats.tsv")   , emit: stats
        tuple val(meta), path("*metabins")            , emit: metabins
        tuple val(meta), path("*metabin_taxonomy.tsv"), emit: taxonomies, optional: true
        path 'versions.yml'                          , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        autometa-binning-summary \\
            --dbdir . \\
            --dbtype ncbi \\
            --binning-main $binning_main \\
            --markers $markers \\
            --metagenome $metagenome \\
            --binning-column $binning_column \\
            --output-stats "${prefix}.metabin_stats.tsv" \\
            --output-taxonomy "${prefix}.metabin_taxonomy.tsv" \\
            --output-metabins "${prefix}.metabins"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
