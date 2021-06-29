
include { hasExtension } from '../modules/local/functions'
// Check mandatory parameters
if (!params.input) exit 1, 'Input samplesheet not specified!'

// Validate input parameters
def hybrid = false
if(hasExtension(params.input, "csv")){
    // just check if long reads are provided
    Channel
        .from(file(params.input))
        .splitCsv(header: true)
        .map { row ->
                if (row.size() == 5) {
                    if (row.long_reads) hybrid = true
                } else {
                    log.error "Input samplesheet contains row with ${row.size()} column(s). Expects 5."
                    System.exit(1)
                }
            }
}
//Workflow.validateWorkflowParams(params, log, hybrid)

def modules = params.modules.clone()

if (params.taxonomy_aware) {
    if (params.single_db_dir == null) {
        error """
        When taxonomy_aware is set to true, You must specify a path to a database directory with --single_db_dir
        """
    }
}

include { ANALYZE_KMERS                                         } from './../modules/local/analyze_kmers'               addParams( options: modules['analyze_kmers']                            )  
include { BIN_CONTIGS                                           } from './../subworkflows/local/bin_contigs.nf'         addParams( binning_options: modules['binning_options'], unclustered_recruitment_options: modules['unclustered_recruitment_options'], binning_summary_options: modules['binning_summary_options']   )                                                
include { GET_SOFTWARE_VERSIONS                                 } from '../modules/local/get_software_versions'         addParams( options: [publish_files : ['csv':'']]                        )
include { HMMER_HMMSEARCH                                       } from './../modules/local/hmmer_hmmsearch.nf'          addParams( options: modules['seqkit_split_options']                     )                                                
include { HMMER_HMMSEARCH_FILTER                                } from './../modules/local/hmmer_hmmsearch_filter.nf'   addParams( options: modules['seqkit_split_options']                     )                                                
include { INPUT_READS                                           } from '../subworkflows/local/input_check'
include { INPUT_CONTIGS                                         } from '../subworkflows/local/input_check'
include { LENGTH_FILTER                                         } from './../modules/local/length_filter'               addParams( options: [publish_files : ['*':'']]                          )
include { MERGE_TSV_WITH_HEADERS as MERGE_SPADES_COVERAGE_TSV   } from './../modules/local/merge_tsv.nf'                addParams( options: modules['seqkit_split_options']                     )                                                
include { MERGE_TSV_WITH_HEADERS as MERGE_MARKERS_TSV           } from './../modules/local/merge_tsv.nf'                addParams( options: modules['merge_markers_options']                     )                                                
include { MERGE_TSV_WITH_HEADERS as MERGE_KMERS_EMBEDDED_TSV    } from './../modules/local/merge_tsv.nf'                addParams( options: modules['seqkit_split_options']                     )                                                
include { MERGE_TSV_WITH_HEADERS as MERGE_KMERS_NORMALIZED_TSV  } from './../modules/local/merge_tsv.nf'                addParams( options: modules['seqkit_split_options']                     )                                                
include { PRODIGAL                                              } from './../modules/nf-core/software/prodigal/main'    addParams( options: modules['prodigal_options']                         )
include { SEARCH_TAXONOMY                                       } from './../subworkflows/local/search_taxonomy'        addParams( diamond_blastp_options: modules['diamond_blastp_options']    )
include { SEQKIT_SPLIT                                          } from './../modules/local/seqkit_split.nf'             addParams( options: modules['seqkit_split_options']  )                                                
include { SPADES_KMER_COVERAGE                                  } from './../modules/local/spades_kmer_coverage'        addParams( options: modules['spades_kmer_coverage']                     )
include { TAXON_ASSIGNMENT                                      } from './../subworkflows/local/taxon_assignment'       addParams( options: modules['taxon_assignment']                         )



workflow AUTOMETA {

    ch_software_versions = Channel.empty()

    INPUT_CONTIGS()
    LENGTH_FILTER(INPUT_CONTIGS.out.metagenome)

    // Split contigs FASTA if running in parallel
    if ( params.parallel_split_fasta ) {
        SEQKIT_SPLIT(LENGTH_FILTER.out.fasta)
        fasta_ch = SEQKIT_SPLIT.out.fasta.transpose()    
    }
    else {        
        fasta_ch = LENGTH_FILTER.out.fasta
    }

    // contig k-mer coverage vs. contig read coverage
    // --------------------------------------------------------------------------------    
    SPADES_KMER_COVERAGE(fasta_ch)
    PRODIGAL(fasta_ch, "gbk")

    // --------------------------------------------------------------------------------
    // OPTIONAL: Run Diamond BLASTp and split data by taxonomy
    // --------------------------------------------------------------------------------        
    if (params.taxonomy_aware) {    
         SEARCH_TAXONOMY(orf_prots_ch)
         TAXON_ASSIGNMENT(fasta_ch, SEARCH_TAXONOMY.out.blastp_table)
         TAXON_ASSIGNMENT.out.taxonomy 
           .set{taxonomy_results}
         ANALYZE_KMERS(TAXON_ASSIGNMENT.out.bacteria)
 
      // TODO KMERS(TAXON_ASSIGNMENT.out.archaea) ... for case of performing binning on archaea
    } else {
        ANALYZE_KMERS(fasta_ch)   
        taxonomy_results = Channel.value( false )
    }
   
    // --------------------------------------------------------------------------------
    // Run hmmsearch and look for marker genes in contig orfs
    // --------------------------------------------------------------------------------  
    HMMER_HMMSEARCH(PRODIGAL.out.amino_acid_fasta, "/home/chase/Documents/github/Autometa/autometa/databases/markers/archaea.single_copy.hmm")
    HMMER_HMMSEARCH.out.domtblout
        .join(PRODIGAL.out.amino_acid_fasta)
        .set{bro}
    HMMER_HMMSEARCH_FILTER(bro)


    // Before binning we need to merge back everything that was run in parallel
    if ( params.parallel_split_fasta ) {
        MERGE_SPADES_COVERAGE_TSV(SPADES_KMER_COVERAGE.out.coverages.groupTuple(),"coverage")
    
        MERGE_SPADES_COVERAGE_TSV.out.merged_tsv
            .set{spades_coverage_merged_tsv_ch}

        MERGE_MARKERS_TSV(HMMER_HMMSEARCH_FILTER.out.markers_tsv.groupTuple(),"markers.tsv")
    
        MERGE_MARKERS_TSV.out.merged_tsv
            .set{markers_tsv_merged_tsv_ch}
 
        MERGE_KMERS_EMBEDDED_TSV(ANALYZE_KMERS.out.embedded.groupTuple(),"kmers.embedded.tsv")
    
        MERGE_KMERS_EMBEDDED_TSV.out.merged_tsv
            .set{kmers_embedded_merged_tsv_ch}
             
        MERGE_KMERS_NORMALIZED_TSV(ANALYZE_KMERS.out.normalized.groupTuple(),"kmers.normalized.tsv")
    
        MERGE_KMERS_NORMALIZED_TSV.out.merged_tsv
            .set{kmers_normalized_tsv_ch}
    }
    else {        
        fasta_ch = LENGTH_FILTER.out.fasta
    }

    BIN_CONTIGS(
        LENGTH_FILTER.out.fasta,
        kmers_embedded_merged_tsv_ch,
        kmers_normalized_tsv_ch,
        spades_coverage_merged_tsv_ch,
        LENGTH_FILTER.out.gc_content,
        markers_tsv_merged_tsv_ch,
        taxonomy_results,
        "cluster"
    )

}
