process LENGTH_TABLE {
    tag "${meta.id}"
    label 'process_low'

    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(metagenome)

    output:
        tuple val(meta), path("lengths.tsv"), emit: lengths
        path  'versions.yml'               , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        #!/usr/bin/env python
        from Bio import SeqIO
        import pandas as pd

        seqs = {record.id: len(record.seq) for record in SeqIO.parse(${metagenome}, "fasta")}
        lengths = pd.Series(seqs, name="length")
        lengths.index.name = "contig"
        lengths.to_csv(lengths.tsv, sep="\t", index=True, header=True)

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
        END_VERSIONS
        """
}
