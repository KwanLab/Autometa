/*
MIT License

Copyright (c) Hadrien Gourlé, Daniel Straub, Sabrina Krakau

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// Import generic module functions
include { hasExtension } from './../../modules/local/functions'

params.options = [:]

workflow INPUT_READS {
    main:
    if(hasExtension(params.input, "csv")){
        // extracts read files from samplesheet CSV and distribute into channels
        Channel
            .from(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                    if (row.size() == 5) {
                        def id = row.sample
                        def group = row.group
                        def sr1 = row.short_reads_1 ? file(row.short_reads_1, checkIfExists: true) : false
                        def sr2 = row.short_reads_2 ? file(row.short_reads_2, checkIfExists: true) : false
                        def lr = row.long_reads ? file(row.long_reads, checkIfExists: true) : false
                        // Check if given combination is valid
                        if (!sr1) exit 1, "Invalid input samplesheet: short_reads_1 can not be empty."
                        if (!sr2 && lr) exit 1, "Invalid input samplesheet: invalid combination of single-end short reads and long reads provided! SPAdes does not support single-end data and thus hybrid assembly cannot be performed."
                        if (!sr2 && !params.single_end) exit 1, "Invalid input samplesheet: single-end short reads provided, but command line parameter `--single_end` is false. Note that either only single-end or only paired-end reads must provided."
                        if (sr2 && params.single_end) exit 1, "Invalid input samplesheet: paired-end short reads provided, but command line parameter `--single_end` is true. Note that either only single-end or only paired-end reads must provided."
                        return [ id, group, sr1, sr2, lr ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 5."
                    }
                }
            .set { ch_input_rows }
        // separate short and long reads
        ch_input_rows
            .map { id, group, sr1, sr2, lr ->
                        def meta = [:]
                        meta.id           = id
                        meta.group        = group
                        meta.single_end   = params.single_end
                        if (params.single_end) 
                            return [ meta, [ sr1] ]
                        else 
                            return [ meta, [ sr1, sr2 ] ]
                }
            .set { ch_raw_short_reads }
        ch_input_rows
            .map { id, group, sr1, sr2, lr ->
                        if (lr) {
                            def meta = [:]
                            meta.id           = id
                            meta.group        = group
                            return [ meta, lr ]
                        }
                }
            .set { ch_raw_long_reads }
    } else {
        Channel
            .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
            .map { row ->
                        def meta = [:]
                        meta.id           = row[0]
                        meta.group        = 0
                        meta.single_end   = params.single_end
                        return [ meta, row[1] ]
                }
            .set { ch_raw_short_reads }
        ch_input_rows = Channel.empty()
        ch_raw_long_reads = Channel.empty()
    }

    // Ensure sample IDs are unique
    ch_input_rows
        .map { id, group, sr1, sr2, lr -> id }
        .toList()
        .map { ids -> if( ids.size() != ids.unique().size() ) {exit 1, "ERROR: input samplesheet contains duplicated sample IDs!" } }

    emit:
    raw_short_reads = ch_raw_short_reads
    raw_long_reads  = ch_raw_long_reads
}



workflow INPUT_CONTIGS {
    main:
        Channel
            .fromPath(params.input)
            .ifEmpty { exit 1, "Cannot find contig fasta file(s)." }
            .map { row ->
                        def meta = [:]
                        meta.id           = row.simpleName
                        return [ meta, row ]
                }
            .set { ch_fasta }
    emit:
        metagenome = ch_fasta
}
