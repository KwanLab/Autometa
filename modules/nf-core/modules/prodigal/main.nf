process PRODIGAL {
    tag "Annotating $meta.id"
    label 'process_low'
    publishDir "${params.outdir}/${meta.id}", mode: params.publish_dir_mode


    conda (params.enable_conda ? "bioconda::prodigal=2.6.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/prodigal:2.6.3--h516909a_2"
    } else {
        container "quay.io/biocontainers/prodigal:2.6.3--h516909a_2"
    }

    input:
    tuple val(meta), path(genome)
    val(output_format)

    output:
    tuple val(meta), path("orfs.${output_format}"), emit: gene_annotations
    tuple val(meta), path("orfs.fna"), emit: nucleotide_fasta
    tuple val(meta), path("orfs.faa"), emit: amino_acid_fasta
    tuple val(meta), path("orfs_all.txt"), emit: all_gene_annotations
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    prodigal -i ${genome} \\
        $options.args \\
        -f $output_format \\
        -d "orfs.fna" \\
        -o "orfs.${output_format}" \\
        -a "orfs.faa" \\
        -s "orfs_all.txt" 

    echo \$(prodigal -v 2>&1) | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p' > ${software}.version.txt
    """
}
