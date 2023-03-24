include { COUNT_KMERS as COUNT            } from '../../modules/local/count_kmers'
include { NORMALIZE_KMERS as NORMALIZE    } from '../../modules/local/normalize_kmers'
include { EMBED_KMERS as EMBED            } from '../../modules/local/embed_kmers'


workflow KMERS {
    take:
        fasta
    main:
        ch_versions = Channel.empty()

        COUNT(fasta)
        ch_versions = ch_versions.mix(COUNT.out.versions)

        NORMALIZE(COUNT.out.counts)
        ch_versions = ch_versions.mix(NORMALIZE.out.versions)

        EMBED(NORMALIZE.out.normalized)
        ch_versions = ch_versions.mix(EMBED.out.versions)

    emit:
        counts = COUNT.out.counts
        normalized = NORMALIZE.out.normalized
        embedded = EMBED.out.embedded
        versions = ch_versions

}

/*
---------------: Test Command :-----------------
nextflow run -resume $HOME/Autometa/subworkflows/local/kmers.nf \\
    --input /path/to/metagenome.fna

*/

workflow {
    ch_fasta = Channel
            .fromPath(params.input)
            .map { row ->
                    def meta = [:]
                    meta.id = row.simpleName
                    return [ meta, row ]
                }

    KMERS(
        ch_fasta
    )
}
