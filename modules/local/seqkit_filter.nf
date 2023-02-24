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
        tuple val(meta), path("*filtered.fna")  , emit: fasta
        tuple val(meta), path("*gc_content.tsv"), emit: gc_content
        path  'versions.yml'                             , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        def metagenomecmd = metagenome.getExtension() == 'gz' ? "gunzip -c $metagenome" : "cat $metagenome"
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        # filter contigs by specified length
        # `seqkit seq -i` "print ID instead of full head"
        ${metagenomecmd} | \\
            seqkit seq -i -j ${task.cpus} -m ${params.length_cutoff} | \\
            seqkit sort -n > "${prefix}.filtered.fna"

        # calculate gc content
        seqkit fx2tab -j ${task.cpus} -n -lg "${prefix}.filtered.fna" > temp

        # Extract columns, create tsv
        awk '{FS="\\t"; OFS="\\t"; print \$1,\$3,\$2}' temp > temp2
        echo -e "contig\\tgc_content\\tlength" | cat - temp2 > "${prefix}.gc_content.tsv"

        # Remove temporary files
        rm temp
        rm temp2


        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            seqkit: \$( seqkit | sed '3!d; s/Version: //' )
        END_VERSIONS
        """
}
