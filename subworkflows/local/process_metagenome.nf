include { CREATE_MOCK                 } from './mock_data'
include { INPUT_CHECK                 } from './input_check'
include { SEQKIT_FILTER               } from '../../modules/local/seqkit_filter'

workflow PROCESS_METAGENOME {

    main:

    ch_versions = Channel.empty()

    assembly_to_locus = Channel.empty()
    assembly_report   = Channel.empty()

    // Set the metagenome and coverage channels
    if (workflow.stubRun){

        CREATE_MOCK()
        ch_versions = ch_versions.mix(CREATE_MOCK.out.versions)

        CREATE_MOCK.out.fasta
            .set{metagenome_ch}

        Channel
            .empty()
            .set{user_provided_coverage_table}

        CREATE_MOCK.out.reads
            .set{reads_ch}

        assembly_to_locus = CREATE_MOCK.out.assembly_to_locus
        assembly_report = CREATE_MOCK.out.assembly_report

    } else {

        Channel
            .fromPath(params.input)
            .set{samplesheet_ch}

        INPUT_CHECK(samplesheet_ch)
        ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

        INPUT_CHECK.out.metagenome
            .set{metagenome_ch}

        INPUT_CHECK.out.coverage
            .set{user_provided_coverage_table}

        INPUT_CHECK.out.reads
            .set{reads_ch}


    }

    SEQKIT_FILTER(
        metagenome_ch
    )
    ch_versions = ch_versions.mix(SEQKIT_FILTER.out.versions)

    SEQKIT_FILTER.out.fasta
        .join(reads_ch)
        .set{combined_contigs_reads}

    emit:
        raw_metagenome_fasta                = metagenome_ch
        filtered_metagenome_fasta           = SEQKIT_FILTER.out.fasta
        user_provided_coverage_table        = user_provided_coverage_table
        reads                               = reads_ch
        filtered_metagenome_fasta_and_reads = combined_contigs_reads
        filtered_metagenome_gc_content      = SEQKIT_FILTER.out.gc_content
        assembly_to_locus                   = assembly_to_locus
        assembly_report                     = assembly_report
        versions                            = ch_versions
}
