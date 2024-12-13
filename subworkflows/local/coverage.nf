include { CALCULATE_COVERAGE   } from './calculate_coverage'
include { SPADES_KMER_COVERAGE } from '../../modules/local/spades_kmer_coverage'

workflow COVERAGE {
    take:
        filtered_metagenome_fasta
        filtered_metagenome_fasta_and_reads
        user_provided_coverage_table

    main:
        // meta.cov_from_assembly.equals('0')

        ch_versions = Channel.empty()

        CALCULATE_COVERAGE(filtered_metagenome_fasta_and_reads)
        ch_versions = ch_versions.mix(CALCULATE_COVERAGE.out.versions)

        SPADES_KMER_COVERAGE (
            filtered_metagenome_fasta,
        )
        ch_versions = ch_versions.mix(SPADES_KMER_COVERAGE.out.versions)

        // https://nextflow-io.github.io/patterns/conditional-process/
        // basically "use input-table coverage, extracted spades coverage, or calculated coverage"
        // TODO: this seems
        user_provided_coverage_table
            .mix(CALCULATE_COVERAGE.out.coverage)
            .mix(SPADES_KMER_COVERAGE.out.coverage)
            .set{coverage_ch}


    emit:
        coverage_ch = coverage_ch
        versions = ch_versions
}
