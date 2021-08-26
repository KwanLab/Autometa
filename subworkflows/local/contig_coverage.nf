params.rev_reads = null
params.fwd_reads = null

params.length_table_options       = [:]
params.align_reads_options        = [:]
params.samtools_viewsort_options  = [:]
params.genome_coverage_options    = [:]

include { LENGTH_TABLE              } from './../../modules/local/length_table.nf'           addParams( options: params.length_table_options                         )
include { ALIGN_READS               } from './../../modules/local/align_reads.nf'            addParams( options: params.align_reads_options                          )
include { SAMTOOLS_VIEW_AND_SORT    } from './../../modules/local/samtools_view_sort.nf'     addParams( samtools_viewsort_options: params.samtools_viewsort_options  )
include { GENOMECOV                 } from './../../subworkflows/local/genome_coverage.nf'   addParams( options: params.genome_coverage_options                      )

workflow CONTIG_COVERAGE {
    take:
        metagenome
        fwd_reads
        rev_reads
        se_reads

    main:
        LENGTH_TABLE(
            metagenome
        )
        ALIGN_READS(
            metagenome,
            fwd_reads,
            rev_reads,
            se_reads
            )
        SAMTOOLS_VIEW_AND_SORT(
            ALIGN_READS.out
            )
        GENOMECOV(
            SAMTOOLS_VIEW_AND_SORT.out,
            LENGTH_TABLE.out
            )

    emit:
        bed = GENOMECOV.out.bed
        coverage = GENOMECOV.out.coverage
}
