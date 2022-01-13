// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LENGTH_TABLE {
    tag "$meta.id"
    label 'process_low'

    publishDir "${meta.id}",
        mode: params.publish_dir_mode,
        saveAs: {
            filename -> saveFiles(
                filename:filename,
                options:params.options,
                publish_dir:getSoftwareName(task.process),
                meta:[:],
                publish_by_meta:[]
            )
        }

    conda (params.enable_conda ? "bioconda::autometa" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("lengths.tsv"), emit: lengths
        path  '*.version.txt'               , emit: version

    script:
        def software = getSoftwareName(task.process)
        """
        #!/usr/bin/env python
        from Bio import SeqIO
        import pandas as pd

        seqs = {record.id: len(record.seq) for record in SeqIO.parse(${metagenome}, "fasta")}
        lengths = pd.Series(seqs, name="length")
        lengths.index.name = "contig"
        lengths.to_csv(lengths.tsv, sep="\t", index=True, header=True)

        echo "TODO" > ${software}.version.txt
        """
}
