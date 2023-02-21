process LENGTH_TABLE {
    tag "${meta.id}"
    label 'process_low'

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

    when:
        task.ext.when == null || task.ext.when

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

        autometa --version | sed -e "s/autometa: //g" > ${software}.version.txt
        """
}
