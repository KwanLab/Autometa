process PRODIGAL {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
    tuple val(meta), path(genome)
    val(output_format)

    output:
    tuple val(meta), path("${prefix}.${output_format}"),    emit: gene_annotations
    tuple val(meta), path("${prefix}.fna"),                 emit: nucleotide_fasta
    tuple val(meta), path("${prefix}.faa"),                 emit: amino_acid_fasta
    tuple val(meta), path("${prefix}_all.txt"),             emit: all_gene_annotations
    path "versions.yml",                                    emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    gzip -cdf ${genome} | prodigal \\
        $args \\
        -f $output_format \\
        -d "${prefix}.fna" \\
        -o "${prefix}.${output_format}" \\
        -a "${prefix}.faa" \\
        -s "${prefix}_all.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
    END_VERSIONS
    """
}
