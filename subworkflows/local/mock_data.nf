include { GET_GENOMES_FOR_MOCK  } from './../../modules/local/get_genomes_for_mock.nf'

process SAMTOOLS_WGSIM {
    // This process is used to create simulated reads from an input FASTA file
    label 'process_low'

    conda "bioconda::samtools=1.13"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0"
    } else {
        container "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }

    input:
    tuple val(meta), path(metagenome)

    output:
    tuple val(meta), path("reads_1.fastq"), path("reads_2.fastq"), emit: reads
    path "versions.yml" , emit: versions

    """
    # https://sarahpenir.github.io/bioinformatics/Simulating-Sequence-Reads-with-wgsim/
    wgsim -1 300 -2 300 -r 0 -R 0 -X 0 -e 0 -N 1 -S 42 ${metagenome} reads_1.fastq reads_2.fastq


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}

workflow CREATE_MOCK {

    main:
        ch_versions = Channel.empty()

        // Download and format fasta files from specfied whole genome assemblies (genomes set from "get_genomes_for_mock" parameter in ~Autometa/conf/modules.config)
        GET_GENOMES_FOR_MOCK()
        ch_versions = ch_versions.mix(GET_GENOMES_FOR_MOCK.out.versions)
        // Format everything with a meta map for use in the main Autometa pipeline
        // see "create_" functions in ~subworkflows/local/input_check.nf
        GET_GENOMES_FOR_MOCK.out.assembly_to_locus
        .map { row ->
                    def meta = [:]
                    meta.id = "mock_data"
                    return [ meta, row ]
            }
        .set { assembly_to_locus }
        GET_GENOMES_FOR_MOCK.out.assembly_report
        .map { row ->
                    def meta = [:]
                    meta.id = "mock_data"
                    return [ meta, row ]
            }
        .set { assembly_report }

        GET_GENOMES_FOR_MOCK.out.metagenome
        .map { row ->
                    def meta = [:]
                    meta.id = "mock_data"
                    return [ meta, row ]
            }
        .set { metagenome }

        // Create fake reads from input genome sequences
        SAMTOOLS_WGSIM(metagenome)
        ch_versions = ch_versions.mix(SAMTOOLS_WGSIM.out.versions)

    emit:
        reads = SAMTOOLS_WGSIM.out.reads
        fasta = metagenome
        assembly_to_locus = assembly_to_locus
        assembly_report = assembly_report
        versions = ch_versions
}
