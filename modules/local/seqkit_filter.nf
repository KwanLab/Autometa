// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SEQKIT_FILTER {
    tag "Removing contigs < ${params.length_cutoff} bp, from ${meta.id}"
    label 'process_high'

    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::seqkit=0.16.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.16.1--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/seqkit:0.16.1--h9ee0642_0"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("${meta.id}.filtered.fna")  , emit: fasta
        tuple val(meta), path("${meta.id}.gc_content.tsv"), emit: gc_content
        path  '*.version.txt'                             , emit: version

    script:
        def software = getSoftwareName(task.process)
        def metagenomecmd = metagenome.getExtension() == 'gz' ? "gunzip -c $metagenome" : "cat $metagenome"
        """
        # filter contigs by specified length
        ${metagenomecmd} | \\
            seqkit seq -j ${task.cpus} -m ${params.length_cutoff} | \\
            seqkit sort -n > "${meta.id}.filtered.fna"

        # calcualte gc content
        seqkit fx2tab -j ${task.cpus} -n -lg "${meta.id}.filtered.fna" > temp

        # Extract columns, create tsv
        awk '{FS="\\t"; OFS="\\t"; print \$1,\$3,\$2}' temp > temp2
        echo -e "contig\\tgc_content\\tlength" | cat - temp2 > "${meta.id}.gc_content.tsv"

        # Remove temporary files
        rm temp
        rm temp2

        seqkit version | sed 's/seqkit v//g' > ${software}.version.txt
        """
}
