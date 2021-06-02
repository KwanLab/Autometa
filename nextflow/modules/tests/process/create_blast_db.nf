#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process CREATE_DIAMOND_DB {
    publishDir params.outdir, pattern: 'diamond_database.dmnd', mode: 'copy'
     input:
        path fasta
    output:
        path 'diamond_database.dmnd'

    """
    diamond makedb \
    --in "${fasta}"\
    --db diamond_database
    """

}
