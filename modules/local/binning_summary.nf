// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

params.taxdump_tar_gz_dir =  [:]

process BINNING_SUMMARY {
    tag "Gathering binning summary for ${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode

    conda (params.enable_conda ? "bioconda::autometa" : null)
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

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-binning-summary \\
            --ncbi $ncbi \\
            --binning-main $binning_main \\
            --markers $markers \\
            --metagenome $metagenome \\
            --binning-column $binning_column \\
            --output-stats "metabin_stats.tsv" \\
            --output-taxonomy "metabin_taxonomy.tsv" \\
            --output-metabins "metabins"

        echo "TODO" > autometa.version.txt
        """
}
