process SPLIT_KINGDOMS {
    tag "Splitting votes into kingdoms for ${meta.id}"
    label 'process_medium'

    conda "bioconda::autometa>=${params.autometa_image_tag}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(assembly), path(votes)
        path taxdump_files // instead of passing to --dbdir, stage and pass '.'

    output:
        tuple val(meta), path("*.tsv.gz")   , emit: contig_taxonomy_tsv
        tuple val(meta), path("*.fna.gz")   , emit: fasta, optional: true
        path  'versions.yml'                , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        autometa-taxonomy \\
            --votes "${votes}" \\
            --output . \\
            --split-rank-and-write superkingdom \\
            --assembly "${assembly}" \\
            --dbdir . \\
            --dbtype ncbi

        # prepend the sample_id to the output files
        for f in *.fna; do mv "\$f" "${meta.id}.\$f"; done

        # rename
        mv taxonomy.tsv "${meta.id}.taxonomy.tsv"

        # compress output files
        gzip -6  "${meta.id}.taxonomy.tsv"
        find ./ -name "*.fna" -type f -exec gzip -6 {} +

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}


