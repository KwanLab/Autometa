#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process GET_ASSEMBLY_SUMMARY {

    output:
      path "assembly_summary_refseq.txt"

    """
    curl -s https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt > assembly_summary_refseq.txt
    """
}


process GET_FTP_DIRS {

    input:
      path assembly_summary_refseq
      val x

    output:
      path "outfile"
    """
    cat "${assembly_summary_refseq}" |\
        grep "${x}" |\
        awk -F '\\t' '{print \$20}' > outfile
    """

}


process DOWNLOAD_MOCK_DATA {

    input:
      path url

    output:
        path "**_genomic.fna.gz", emit: nucleotide
        path "**_protein.faa.gz", emit: protein

    """
    cat outfile | sed 's,ftp://,rsync://,' | xargs -n 1 -I {} rsync -am --exclude='*_rna_from_genomic.fna.gz' --exclude='*_cds_from_genomic.fna.gz' --include="*_genomic.fna.gz" --include="*_protein.faa.gz" --include='*/' --exclude='*' {} .
    """
}

process WRITE_FILE_APPEND_COV {
    publishDir params.outdir, mode: 'copy'

    input:
        path x
        val y

    output:
        path "${y}"  , emit: fasta

    """
    cat "${x}" | awk '/^>/ {\$0=\$1} 1' | sed 's/>.*/&_length_1_cov_1/' > "${y}"
    """



}


assemblies = Channel.fromList(
    [
        "GCF_008124965.1"
    ]
)

workflow CREATE_MOCK {

    main:

        GET_ASSEMBLY_SUMMARY()
        GET_FTP_DIRS(
            GET_ASSEMBLY_SUMMARY.out,
            assemblies.flatten()
            )
        DOWNLOAD_MOCK_DATA(GET_FTP_DIRS.out)
        WRITE_FILE_APPEND_COV(
            DOWNLOAD_MOCK_DATA.out.nucleotide.splitFasta(by:1).collectFile(),
            "mock_metagenome.fna"
        )

        WRITE_FILE_APPEND_COV.out.fasta
        .map { row ->
                    def meta = [:]
                    meta.id           = row.simpleName
                    return [ meta, row ]
            }
        .set { ch_fasta }

    emit:
        fasta = ch_fasta


}

