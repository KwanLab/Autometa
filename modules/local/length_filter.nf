// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process LENGTH_FILTER {
    tag "Removing ${meta.id} contigs less than ${params.length_cutoff}"
    label 'process_high'
   
    publishDir "${params.interim_dir_internal}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
    
    
    conda (params.enable_conda ? "bioconda::seqkit=0.15.0" : null)

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/seqkit:0.15.0--0"
    } else {
        container "quay.io/biocontainers/seqkit:0.15.0--0"
    }

    input:
        tuple val(meta), path(metagenome)
  
    output:
      tuple val(meta), path("${meta.id}.filtered.fna"), emit: fasta
      tuple val(meta), path("${meta.id}.gc_content.tsv"), emit: gc_content
    //  tuple val(meta), path("${meta.id}.stats.tsv"), emit: stats
      
  
    script:
    def software = getSoftwareName(task.process)

   """
    zcat ${metagenome} |\
        seqkit seq -j ${task.cpus} -m ${params.length_cutoff} > "${meta.id}.filtered.fna"

    cat "${meta.id}.filtered.fna" |\
        seqkit fx2tab -j ${task.cpus} -n -lg > temp
    
    awk '{FS="\\t"; OFS="\\t"; print \$1,\$3,\$2}' temp > temp2
    echo -e "contig\\tgc_content\\tlength" | cat - temp2 > "${meta.id}.gc_content.tsv"
    rm temp
    rm temp2
    """

    // python code fails if input file ends with .filtered.fna 
    //if ( "${meta.id}" == "${meta.id}.filtered.fna" )
    //"""
    //autometa-length-filter \
    //  --assembly $metagenome \
    //  --cutoff ${params.length_cutoff} \
    //  --output-fasta "${meta.id}.autometa_filtered.fna" \
    //  --output-stats ${meta.id}.stats.tsv \
    //  --output-gc-content ${meta.id}.gc_content.tsv
    //"""
    //else 
    //"""
    //autometa-length-filter \
    //  --assembly $metagenome \
    //  --cutoff ${params.length_cutoff} \
    //  --output-fasta "${meta.id}.filtered.fna" \
    //  --output-stats ${meta.id}.stats.tsv \
    //  --output-gc-content ${meta.id}.gc_content.tsv
  //
    //"""
}
