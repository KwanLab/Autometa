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




process WRITE_FILE {
    publishDir params.outdir, mode: 'copy'

    input:
        path x
        val y
    
    output:
        path "${y}"

    """
    cat "${x}" > "${y}"
    """



}
