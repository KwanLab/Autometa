
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




include { GET_SOFTWARE_VERSIONS                               } from '../modules/local/get_software_versions'       addParams( options: [publish_files : ['csv':'']] )

include { INPUT_READS         } from '../subworkflows/local/input_check'
include { INPUT_CONTIGS         } from '../subworkflows/local/input_check'

include { LENGTH_FILTER } from './../modules/local/length_filter' addParams( options: [publish_files : ['*':'']] )
include { PRODIGAL } from './../modules/nf-core/software/prodigal/main'


//
include { ANALYZE_KMERS } from './../modules/local/analyze_kmers' addParams( options: modules['analyze_kmers'] )
include { TAXON_ASSIGNMENT } from './../subworkflows/local/taxon_assignment' addParams( options: modules['taxon_assignment'] )
include { SEARCH_TAXONOMY } from './../subworkflows/local/search_taxonomy' addParams( options: modules['search_taxonomy'] )

include { SPADES_KMER_COVERAGE } from './../modules/local/spades_kmer_coverage' addParams( options: modules['spades_kmer_coverage'] )

include { MARKERS } from './../modules/local/markers' addParams( options: modules['markers'] )

include { BIN_CONTIGS } from './../subworkflows/local/bin_contigs.nf' addParams( binning_options: modules['binning_options'], unclustered_recruitment_options: modules['unclustered_recruitment_options'], binning_summary_options: modules['binning_summary_options']   )                                                



// if (params.gtdb) {
//     Channel
//         .value(file( "${params.gtdb}" ))
//         .set { ch_gtdb }
// } else {
//     ch_gtdb = Channel.empty()
// }


//def busco_failed_bins = [:]


metagenome_ch = ''


workflow AUTOMETA {

    ch_software_versions = Channel.empty()

    INPUT_CONTIGS()
    LENGTH_FILTER(INPUT_CONTIGS.out.metagenome)

    // contig k-mer coverage vs. contig read coverage
    // --------------------------------------------------------------------------------    
    SPADES_KMER_COVERAGE(LENGTH_FILTER.out.fasta)
    

    if ( params.parallel_split_fasta ) {

           LENGTH_FILTER.out.fasta.splitFasta(file:"${params.store_split_fasta_in_ram}", size:4.MB) //TODO: Add parameter for number of splits
             .set{filtered_ch}      

           PRODIGAL(filtered_ch, "gbk")

           PRODIGAL.out.amino_acid_fasta
             .set{orf_prots_parallel_ch}

           PRODIGAL.out.amino_acid_fasta.collectFile(cache:false)
             .set{orf_prots_ch}

    } else {
       
           PRODIGAL(LENGTH_FILTER.out.fasta, "gbk")

           PRODIGAL.out.amino_acid_fasta
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

      MARKERS(
        orf_prots_parallel_ch
        )    
      
      MARKERS.out.groupTuple().view()
      .collectFile(
        cache:true,
        keepHeader:true
        )
        .view()
        .map { meta, reads -> [ meta.id, meta, reads ] }
        .set{
          markers_out_ch
          }

    }
    else {
      MARKERS(orf_prots_ch)
      MARKERS.out.set{markers_out_ch}
  
    }


    // --------------------------------------------------------------------------------
    // Bin contigs
    // --------------------------------------------------------------------------------  

 

    BIN_CONTIGS(
      LENGTH_FILTER.out.fasta,
      ANALYZE_KMERS.out.embedded,
      ANALYZE_KMERS.out.normalized,
      SPADES_KMER_COVERAGE.out,
      LENGTH_FILTER.out.gc_content,
      markers_out_ch,
      taxonomy_results,
      "cluster"
    )

}
