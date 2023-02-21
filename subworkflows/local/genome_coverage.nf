include { BEDTOOLS_GENOMECOV } from './../../modules/nf-core/modules/bedtools/genomecov.nf'

workflow GENOME_COVERAGE {
    take:
        bam         // channel: [ val(meta), path(bam) ]
        lengths     // channel: [ val(meta), path(lengths) ]  // https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format

    main:
        bedtools_input_ch = bam.combine(lengths)

        BEDTOOLS_GENOMECOV (
            bedtools_input_ch
            )

        bam.out.bed
            .combine(lengths)
            .combine(BEDTOOLS_GENOMECOV.out.bed)
            .set{parse_bed_input_ch}

        PARSE_BED (
            parse_bed_input_ch
        )

    emit:
        bed = BEDTOOLS_GENOMECOV.out.bed
        coverage = PARSE_BED.out.coverage
}
