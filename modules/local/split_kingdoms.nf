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
        val dbtype

    output:
        tuple val(meta), path("${dbtype}/*.taxonomy.tsv")      , emit: taxonomy
        tuple val(meta), path("${dbtype}/*.fna")               , emit: fna
        tuple val(meta), path("${dbtype}/*.unclassified.fna")  , emit: unclassified_fna, optional: true
        path  'versions.yml'                                   , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        mkdir ${dbtype}
        autometa-taxonomy \\
            --votes "${votes}" \\
            --output "./${dbtype}" \\
            --split-rank-and-write superkingdom \\
            --assembly "${assembly}" \\
            --dbdir . \\
            --dbtype ${dbtype}

        # prefix all files in temp with the prefix
        for file in ${dbtype}/*; do
            mv "\$file" "${dbtype}/${prefix}.\$(basename \$file)"
        done

        # Move .unclassified.fna files to a separate location for separate emitting
        mkdir -p ${dbtype}_unclassified_fna

        for file in ${dbtype}/${prefix}.unclassified.*; do
            if [ -e "\$file" ]; then
                mv "\$file" ${dbtype}_unclassified_fna/
            fi
        done

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
