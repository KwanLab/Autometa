include { initOptions; saveFiles; getSoftwareName } from './functions'

params.get_genomes_for_mock = [:]
include { GET_GENOMES_FOR_MOCK  } from './../../modules/local/get_genomes_for_mock.nf' addParams( options: params.get_genomes_for_mock )






process SAMTOOLS_WGSIM {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::samtools=1.13" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/samtools:1.12--hd5e65b6_0"
    } else {
        container "quay.io/biocontainers/samtools:1.12--hd5e65b6_0"
    }

    input:
    path fasta

    output:
    path("*.fastq"), emit: fastq
    path "*.version.txt"          , emit: version


    """
    # https://sarahpenir.github.io/bioinformatics/Simulating-Sequence-Reads-with-wgsim/
    wgsim -1 300 -2 300 -r 0 -R 0 -X 0 -e 0 ${fasta} reads_1.fastq reads_2.fastq

    echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' > samtools.version.txt
    """
}


// "GCF_013307045.1" has no markers
workflow CREATE_MOCK {

    main:
        GET_GENOMES_FOR_MOCK()

        SAMTOOLS_WGSIM(GET_GENOMES_FOR_MOCK.out.metagenome)

        GET_GENOMES_FOR_MOCK.out.fake_spades_coverage
        .map { row ->
                    def meta = [:]
                    meta.id = 'mock_data'
                    return [ meta, row ]
            }
        .set { ch_fasta }
        GET_GENOMES_FOR_MOCK.out.assembly_to_locus
        .map { row ->
                    def meta = [:]
                    meta.id = 'mock_data'
                    return [ meta, row ]
            }
        .set { assembly_to_locus }
        GET_GENOMES_FOR_MOCK.out.assembly_report
        .map { row ->
                    def meta = [:]
                    meta.id = 'mock_data'
                    return [ meta, row ]
            }
        .set { assembly_report }

    emit:
        fasta = ch_fasta
        reads = SAMTOOLS_WGSIM.out.fastq
        assembly_to_locus = assembly_to_locus
        assembly_report =   assembly_report
}
