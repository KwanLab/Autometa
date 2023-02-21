
//TODO: These don't map to anything
params.fwd_reads = null
params.rev_reads = null
params.se_reads = null

include { ALIGN_READS               } from '../../modules/local/align_reads'
include { SAMTOOLS_VIEW_AND_SORT    } from '../../modules/local/samtools_view_sort'
include { BEDTOOLS_GENOMECOV        } from '../../modules/local/bedtools_genomecov'
include { PARSE_BED                 } from '../../modules/local/parse_bed'

workflow CONTIG_COVERAGE {
    take:
        metagenome_reads_ch

    main:
        ALIGN_READS(
            metagenome_reads_ch
        )
        SAMTOOLS_VIEW_AND_SORT(
            ALIGN_READS.out.sam
        )
        BEDTOOLS_GENOMECOV(
            SAMTOOLS_VIEW_AND_SORT.out.bam
        )
        PARSE_BED(BEDTOOLS_GENOMECOV.out.bed)

    emit:
        sam = ALIGN_READS.out.sam
        bam = SAMTOOLS_VIEW_AND_SORT.out.bam
        bed = BEDTOOLS_GENOMECOV.out.bed
        coverage = PARSE_BED.out.coverage
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
