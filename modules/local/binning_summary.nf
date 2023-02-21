process BINNING_SUMMARY {
    tag "Gathering binning summary for ${meta.id}"
    label 'process_high'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(binning_main), path(markers), path(metagenome)
        val(binning_column)
        path(ncbi)

    output:
        tuple val(meta), path("metabin_stats.tsv")   , emit: stats
        tuple val(meta), path("metabins")            , emit: metabins
        tuple val(meta), path("metabin_taxonomy.tsv"), emit: taxonomies, optional: true
        path  '*.version.txt'                        , emit: version

    when:
        task.ext.when == null || task.ext.when

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-binning-summary \\
            --dbdir $ncbi \\
            --dbtype ncbi \\
            --binning-main $binning_main \\
            --markers $markers \\
            --metagenome $metagenome \\
            --binning-column $binning_column \\
            --output-stats "metabin_stats.tsv" \\
            --output-taxonomy "metabin_taxonomy.tsv" \\
            --output-metabins "metabins"

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
