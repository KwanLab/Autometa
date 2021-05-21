#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { LENGTH_FILTER } from './autometa_core/utilities/process/metagenome_length_filter.nf'
include { ANALYZE_KMERS } from './autometa_core/utilities/process/analyze_kmers.nf'
include { KMER_COVERAGE } from './autometa_core/utilities/process/kmer_coverage.nf'
include { PRODIGAL_ORFS } from './autometa_core/utilities/process/orf_calling.nf'
include { MARKERS } from './autometa_core/utilities/process/markers.nf'
include { TAXON_ASSIGNMENT } from './contig_split_by_taxonomy/taxonomy_workflow.nf'
include { SEARCH_TAXONOMY } from './contig_search_taxonomy/search_taxonomy.nf'
include { BIN_CONTIGS } from './autometa_core/contig_binning/binning_workflow.nf'


workflow AUTOMETA {
  take:
    metagenome

  main:
    // Perform various annotations on provided metagenome
    LENGTH_FILTER(metagenome)
    // --------------------------------------------------------------------------------
    // k-mer coverage vs. read coverage
    // --------------------------------------------------------------------------------    
    KMER_COVERAGE(LENGTH_FILTER.out.fasta)
    
    // --------------------------------------------------------------------------------
    // Run Prodigal to obtain open reading frames
    // --------------------------------------------------------------------------------        
    if ( params.parallel_split_fasta ) {
      
      LENGTH_FILTER.out.fasta.splitFasta(file:"${params.store_split_fasta_in_ram}", size:4.MB) //TODO: Add parameter for number of splits
        .set{filtered_ch}      
      
      PRODIGAL_ORFS(filtered_ch)
      
      PRODIGAL_ORFS.out.prots
        .set{orf_prots_parallel_ch}

      PRODIGAL_ORFS.out.prots.collectFile(cache:false)
        .set{orf_prots_ch}

    }
    else {

      PRODIGAL_ORFS(LENGTH_FILTER.out.fasta) 
      PRODIGAL_ORFS.out.prots
        .set{orf_prots_ch}

    }
    // --------------------------------------------------------------------------------
    // OPTIONAL: Run Diamond BLASTp and split data by taxonomy
    // --------------------------------------------------------------------------------        
    if (params.taxonomy_aware) {    
      SEARCH_TAXONOMY(orf_prots_ch)
      TAXON_ASSIGNMENT(LENGTH_FILTER.out.fasta, SEARCH_TAXONOMY.out.blastp_table)
      
      TAXON_ASSIGNMENT.out.taxonomy 
        .set{taxonomy_results}
      ANALYZE_KMERS(TAXON_ASSIGNMENT.out.bacteria)

      // TODO KMERS(TAXON_ASSIGNMENT.out.archaea) ... for case of performing binning on archaea
    } else {
    ANALYZE_KMERS(LENGTH_FILTER.out.fasta)   
    taxonomy_results = Channel.value( false )
    }

    // --------------------------------------------------------------------------------
    // Run hmmscan and look for marker genes in contig orfs
    // --------------------------------------------------------------------------------        

    if ( params.parallel_split_fasta ) {

      MARKERS(orf_prots_parallel_ch)    
      MARKERS.out.collectFile(cache:true, keepHeader:true)
        .set{markers_out_ch}

    }
    else {
    MARKERS(orf_prots_ch)
    MARKERS.out.set{markers_out_ch}
  
    }
    // --------------------------------------------------------------------------------
    // Bin contigs
    // --------------------------------------------------------------------------------  
    BIN_CONTIGS(LENGTH_FILTER.out.fasta, ANALYZE_KMERS.out.embedded, ANALYZE_KMERS.out.normalized, KMER_COVERAGE.out, LENGTH_FILTER.out.gc_content, markers_out_ch, taxonomy_results, "cluster")
}
