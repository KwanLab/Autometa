
process SEQKIT_FILTER {
    tag "Removing contigs < ${params.length_cutoff} bp, from ${meta.id}"
    label 'process_high'

    conda "bioconda::seqkit=0.16.1"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("*filtered.fna.gz")   , emit: fasta
        tuple val(meta), path("*gc_content.tsv.gz") , emit: gc_content
        path  'versions.yml'                        , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        # filter contigs by specified length
        gzip -fdc ${metagenome} | \\
            seqkit seq ${args} -j ${task.cpus} -m ${params.length_cutoff} | \\
            seqkit sort ${args2} -n > "${prefix}.filtered.fna"

        # calculate gc content
        seqkit fx2tab -j ${task.cpus} -n -lg "${prefix}.filtered.fna" > temp

        # Extract columns, create tsv with columns and column order the autometa python code expects
        awk '{FS="\\t"; OFS="\\t"; print \$1,\$3,\$2}' temp > temp2
        echo -e "contig\\tgc_content\\tlength" | cat - temp2 > "${prefix}.gc_content.tsv"

        # Remove temporary files
        rm temp
        rm temp2

        gzip -6  "${prefix}.filtered.fna"
        gzip -6  "${prefix}.gc_content.tsv"

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        END_VERSIONS
        """
}
