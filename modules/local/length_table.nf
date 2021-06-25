// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LENGTH_TABLE {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jason-c-kwan/autometa:nfcore"
    }

    input:
    tuple val(meta), path(metagenome)

    output:   
    tuple val(meta), path("${meta.id}.lengths.tsv"), emit: bam

    script:
    def software = getSoftwareName(task.process)
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    import pandas as pd

    seqs = {record.id: len(record) for record in SeqIO.parse(${metagenome}, "fasta")}
    lengths = pd.Series(seqs, name="length")
    lengths.index.name = "contig"
    lengths.to_csv(${meta.id}.lengths.tsv, sep="\t", index=True, header=True)

    # echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > ${software}.version.txt
    """
}
