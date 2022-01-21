//
// Check input samplesheet and get read channels
//

nextflow.enable.dsl=2

params.options = [:]

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check' addParams( options: params.options )

workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK ( samplesheet )

        // reads channel
        SAMPLESHEET_CHECK.out
            .splitCsv ( header:true, sep:',' )
            .map { create_fastq_channel(it) }
            .set { reads }
        // metagenome channel
        SAMPLESHEET_CHECK.out
            .splitCsv ( header:true, sep:',' )
            .map { create_metagenome_channel(it) }
            .set { metagenome }

        // coverage channel
        SAMPLESHEET_CHECK.out
            .splitCsv ( header:true, sep:',' )
            .map { create_coverage_channel(it) }
            .set { coverage }
    emit:
        reads      // channel: [ val(meta), [ reads ] ]
        metagenome // channel: [ val(meta), [ assembly ]]
        coverage   // channel: [ val(meta), [ coverage ]]
}


// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id               = row.sample
    meta.single_end       = row.single_end.toBoolean()
    meta.cov_from_headers = row.cov_from_contig_headers.toBoolean()

    def array = []
    if (meta.cov_from_headers || file(row.coverage_tab).exists()) {
        return
    }
    if (!meta.cov_from_headers and !file(row.fastq_1).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist and this sample was specified to compute coverage using reads!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        array = [ meta, file(row.fastq_1), "0" ]
    } else {
        if (!file(row.fastq_2).exists() and !meta.cov_from_headers) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist and this sample was specified to compute coverage using reads!\n${row.fastq_2}"
        }
        array = [ meta, file(row.fastq_1), file(row.fastq_2) ]
    }
    return array
}

// Function to get list of [ meta, assembly ]
def create_metagenome_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id               = row.sample
    meta.single_end       = row.single_end.toBoolean()
    meta.cov_from_headers = row.cov_from_contig_headers.toBoolean()

    def array = []
    if (!file(row.assembly).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Metagenome assembly file does not exist!\n${row.assembly}"
    }
    array = [ meta, [ file(row.assembly) ] ]
    return array
}

// Function to get list of [ meta, coverage ]
def create_coverage_channel(LinkedHashMap row) {
    def meta = [:]
    meta.id               = row.sample
    meta.single_end       = row.single_end.toBoolean()
    meta.cov_from_headers = row.cov_from_contig_headers.toBoolean()

    def array = []
    if (!meta.cov_from_headers and !file(row.fastq_1).exists() and !file(row.coverage_tab).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Coverage file does not exist!\n${row.coverage_tab}"
    }
    if (meta.cov_from_headers || !file(row.coverage_tab).exists()) {
        // e.g. get coverage from headers or read alignments...
        // Don't create a coverage tab channel for the sample if it was not provided...
        return
    } else {
        array = [ meta, [ file(row.coverage_tab) ] ]
        return array
    }
}


/*
---------------: Test Command :-----------------
nextflow run $HOME/Autometa/subworkflows/local/input_check.nf \\
    --input test_samplesheet.csv

/home/evan/autometa_test_data/autometa_test_data/simulated_communities/78.125Mbp/

forward_reads.fastq.gz
reverse_reads.fastq.gz
metagenome.fna.gz

*/

workflow {
    main:
        samplesheet_ch = Channel.fromPath(params.input)
        INPUT_CHECK(samplesheet_ch)
    emit:
        reads = INPUT_CHECK.out.reads
        metagenome = INPUT_CHECK.out.metagenome
}
