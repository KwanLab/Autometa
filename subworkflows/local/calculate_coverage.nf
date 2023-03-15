
//TODO: These don't map to anything
params.fwd_reads = null
params.rev_reads = null
params.se_reads = null

include { ALIGN_READS               } from '../../modules/local/align_reads'
include { BEDTOOLS_GENOMECOV        } from '../../modules/local/bedtools_genomecov'
include { PARSE_BED                 } from '../../modules/local/parse_bed'

workflow CALCULATE_COVERAGE {
    take:
        metagenome_reads_ch

    main:
        ch_versions = Channel.empty()

        ALIGN_READS(
            metagenome_reads_ch
        )
        ch_versions = ch_versions.mix(ALIGN_READS.out.versions)

        BEDTOOLS_GENOMECOV(
            ALIGN_READS.out.bam
        )
        ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

        PARSE_BED(BEDTOOLS_GENOMECOV.out.bed)
        ch_versions = ch_versions.mix(PARSE_BED.out.versions)

    emit:
        bam = ALIGN_READS.out.bam
        bed = BEDTOOLS_GENOMECOV.out.bed
        coverage = PARSE_BED.out.coverage
        versions = ch_versions

}

/*
---------------: Test Command :-----------------
nextflow run -resume $HOME/Autometa/subworkflows/local/contig_coverage.nf \\
    --publish_dir_mode copy \\
    --fwd_reads '/path/to/fwd_reads.fastq.gz' \\
    --rev_reads '/path/to/rev_reads.fastq.gz'

*/

workflow {
    coverage_inputs_ch = Channel
            .fromPath(params.input)
            .map { row ->
                    def meta = [:]
                    meta.id = row.simpleName
                    return [ meta, row ]
                }

    CONTIG_COVERAGE(
        coverage_inputs_ch
    )
}
