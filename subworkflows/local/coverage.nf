include { CALCULATE_COVERAGE   } from './calculate_coverage'
include { SPADES_KMER_COVERAGE } from '../../modules/local/spades_kmer_coverage'

workflow COVERAGE {
    take:
        filtered_metagenome_fasta
        filtered_metagenome_fasta_and_reads
        user_provided_coverage_table

    main:

        ch_versions = Channel.empty()


        CALCULATE_COVERAGE(filtered_metagenome_fasta_and_reads)
        ch_versions = ch_versions.mix(CALCULATE_COVERAGE.out.versions)

        filtered_metagenome_fasta
            .filter( {meta, file -> meta.cov_from_assembly == 'spades' })
            .set {ch_for_spades_cov}

        SPADES_KMER_COVERAGE (
            ch_for_spades_cov,
        )
        ch_versions = ch_versions.mix(SPADES_KMER_COVERAGE.out.versions)

        // https://nextflow-io.github.io/patterns/conditional-process/
        // basically "use input-table coverage, extracted spades coverage, or calculated coverage"
        // TODO: this seems like it should choose one or the other, not mix?
        user_provided_coverage_table
            .mix(CALCULATE_COVERAGE.out.coverage)
            .mix(SPADES_KMER_COVERAGE.out.coverage)
            .set{ch_coverage}


    emit:
        ch_coverage = ch_coverage
        versions = ch_versions
}
