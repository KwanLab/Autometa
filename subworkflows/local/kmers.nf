include { COUNT_KMERS as COUNT            } from '../../modules/local/count_kmers'     addParams( options: params.count_kmers_options      )
include { NORMALIZE_KMERS as NORMALIZE    } from '../../modules/local/normalize_kmers' addParams( options: params.normalize_kmers_options  )
include { EMBED_KMERS as EMBED            } from '../../modules/local/embed_kmers'     addParams( options: params.embed_kmers_options      )


workflow KMERS {
    take:
        fasta
    main:
        COUNT(fasta)
        NORMALIZE(COUNT.out.counts)
        EMBED(NORMALIZE.out.normalized)
    emit:
        counts = COUNT.out.counts
        normalized = NORMALIZE.out.normalized
        embedded = EMBED.out.embedded
}

/*
---------------: Test Command :-----------------
nextflow run -resume $HOME/Autometa/subworkflows/local/kmers.nf \\
    --input /path/to/metagenome.fna

*/

workflow {
    fasta_ch = Channel
            .fromPath(params.input)
            .map { row ->
                    def meta = [:]
                    meta.id = row.simpleName
                    return [ meta, row ]
                }

    KMERS(
        fasta_ch
    )
}
