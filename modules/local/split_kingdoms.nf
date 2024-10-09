process SPLIT_KINGDOMS {
    tag "Splitting votes into kingdoms for ${meta.id}"
    label 'process_medium'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(assembly), path(votes)
        path taxdump_files // instead of passing to --dbdir, stage and pass '.'

    output:
        tuple val(meta), path("*.taxonomy.tsv"), emit: taxonomy
        tuple val(meta), path("*.bacteria.fna"), emit: bacteria, optional: true
        tuple val(meta), path("*.archaea.fna") , emit: archaea,  optional: true
        tuple val(meta), path("*.fna")                 , emit: kingdoms, optional: true
        path  'versions.yml'                           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir temp
        autometa-taxonomy \\
            --votes "${votes}" \\
            --output "./temp" \\
            --split-rank-and-write superkingdom \\
            --assembly "${assembly}" \\
            --dbdir . \\
            --dbtype ncbi

        # prefix all files in temp with the prefix
        for file in temp/*; do
            mv "\$file" "${prefix}.\$(basename \$file)"
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
