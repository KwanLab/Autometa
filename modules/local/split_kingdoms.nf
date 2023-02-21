process SPLIT_KINGDOMS {
    tag "Splitting votes into kingdoms for ${meta.id}"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(assembly), path(votes)
        path(ncbi_tax_dir)

    output:
        tuple val(meta), path("taxonomy.tsv"), emit: taxonomy
        tuple val(meta), path("bacteria.fna"), emit: bacteria, optional: true
        tuple val(meta), path("archaea.fna") , emit: archaea,  optional: true
        tuple val(meta), path("*.fna")       , emit: kingdoms, optional: true
        path  '*.version.txt'                , emit: version

    when:
        task.ext.when == null || task.ext.when

    script:
        def software = getSoftwareName(task.process)
        """
        autometa-taxonomy \\
            --votes "${votes}" \\
            --output . \\
            --split-rank-and-write superkingdom \\
            --assembly "${assembly}" \\
            --dbdir "${ncbi_tax_dir}" \\
            --dbtype ncbi

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
