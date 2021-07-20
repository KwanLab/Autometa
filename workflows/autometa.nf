/*
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
 */
//Workflow.validateWorkflowParams(params, log, hybrid)

def modules = params.modules.clone()


def check_for_file(path) {
    return
}


// check if user wants to separate contigs based on taxonomy before binning

    if (params.single_db_dir) {
        nr_dmnd_dir = "${params.single_db_dir}"
        prot_accession2taxid_gz_dir = "${params.single_db_dir}"
        taxdump_tar_gz_dir = "${params.single_db_dir}"
    }
    // instead of doing '(!params.single_db_dir && params.nr_dmnd_dir)'
    // just override e.g. 'nr_dmnd_location' so users can set
    // 'single_db_dir' but also set individual other db paths if they have them
    // e.g. if they have nr.dmnd but not the other files.
    if (!params.single_db_dir && params.nr_dmnd_dir) {
        nr_dmnd_dir = "${params.nr_dmnd_dir}"
    }
    if (!params.single_db_dir && params.prot_accession2taxid_gz_dir) {
        prot_accession2taxid_gz_dir = "${params.taxdump_tar_gz_dir}" // currently the python needs it to be here
    }
    if (!params.single_db_dir && params.taxdump_tar_gz_dir) {
        taxdump_tar_gz_dir = "${params.taxdump_tar_gz_dir}"
    }
    if (params.large_downloads_permission) {
        // TODO: check if files already exist, if they do fail the pipeline early at this stage
    } else {
        // TODO: check if files exist, if they don't fail the pipeline early at this stage
    }


include { ANALYZE_KMERS                                         } from './../modules/local/analyze_kmers'               addParams( options: modules['analyze_kmers']                            )
include { BIN_CONTIGS                                           } from './../subworkflows/local/bin_contigs.nf'         addParams( binning_options: modules['binning_options'], unclustered_recruitment_options: modules['unclustered_recruitment_options'], binning_summary_options: modules['binning_summary_options'], prot_accession2taxid_gz_dir: prot_accession2taxid_gz_dir )
include { GET_SOFTWARE_VERSIONS                                 } from '../modules/local/get_software_versions'         addParams( options: [publish_files : ['csv':'']]                        )
include { HMMER_HMMSEARCH                                       } from './../modules/local/hmmer_hmmsearch.nf'          addParams( options: modules['hmmsearch_options']                        )
include { HMMER_HMMSEARCH_FILTER                                } from './../modules/local/hmmer_hmmsearch_filter.nf'   addParams( options: modules['hmmsearch_filter_options']                 )
include { INPUT_CONTIGS                                         } from '../subworkflows/local/input_check'              addParams(                                                              )
include { CREATE_MOCK                                           } from '../subworkflows/local/mock_data.nf'             addParams()
include { SEQKIT_FILTER                                         } from './../modules/local/seqkit_filter'               addParams( options: [publish_files : ['*':'']]                          )
include { MERGE_TSV_WITH_HEADERS as MERGE_SPADES_COVERAGE_TSV   } from './../modules/local/merge_tsv.nf'                addParams( options: modules['spades_kmer_coverage']                     )
include { MERGE_TSV_WITH_HEADERS as MERGE_HMMSEARCH             } from './../modules/local/merge_tsv.nf'                addParams( options: modules['merge_hmmsearch_options']                  )
include { PRODIGAL                                              } from './../modules/nf-core/modules/prodigal/main'     addParams( options: modules['prodigal_options']                         )
include { SEQKIT_SPLIT                                          } from './../modules/local/seqkit_split.nf'             addParams( options: modules['seqkit_split_options'], num_splits: params.num_splits                     )
include { SPADES_KMER_COVERAGE                                  } from './../modules/local/spades_kmer_coverage'        addParams( options: modules['spades_kmer_coverage']                     )
include { TAXON_ASSIGNMENT                                      } from './../subworkflows/local/taxon_assignment'       addParams( options: modules['taxon_assignment'], majority_vote_options: modules['majority_vote_options'], split_kingdoms_options: modules['split_kingdoms_options']                         )
include { PREPARE_NR_DB                                         } from './../subworkflows/local/prepare_nr.nf'          addParams( debug:params.debug, diamond_makedb_options: modules['diamond_makedb_options'], nr_dmnd_dir: nr_dmnd_dir)
include { PREPARE_TAXONOMY_DATABASES                            } from './../subworkflows/local/prepare_ncbi_taxinfo.nf' addParams( debug:params.debug, taxdump_tar_gz_dir: taxdump_tar_gz_dir, prot_accession2taxid_gz_dir: prot_accession2taxid_gz_dir)
include { DIAMOND_BLASTP                                        } from './../modules/local/diamond_blastp.nf'           addParams( options: modules['diamond_blastp_options']                   )
include { MERGE_FASTA as MERGE_PRODIGAL                         } from './../modules/local/merge_fasta.nf'              addParams()
include { MARKERS                                               } from './../modules/local/markers.nf'                  addParams( options: modules['seqkit_split_options']                     )

workflow AUTOMETA {
    ch_software_versions = Channel.empty()

if (params.mock_test){

    CREATE_MOCK()

    CREATE_MOCK.out.fasta
        .set{input_ch}

} else {

    INPUT_CONTIGS()

    INPUT_CONTIGS.out.metagenome
        .set{input_ch}
}

    SEQKIT_FILTER(input_ch)

    // Split contigs FASTA if running in parallel
    if ( params.parallel_split_fasta ) {
        SEQKIT_SPLIT(SEQKIT_FILTER.out.fasta)
        fasta_ch = SEQKIT_SPLIT.out.fasta.transpose()
    }
    else {
        fasta_ch = SEQKIT_FILTER.out.fasta
    }

/*
 * -------------------------------------------------
 *  Find coverage, currently only pulling from SPADES output
 * -------------------------------------------------
 */
    SPADES_KMER_COVERAGE(fasta_ch)


/*
 * -------------------------------------------------
 *  Find open reading frames with Prodigal
 * -------------------------------------------------
 */
    PRODIGAL(fasta_ch, "gbk")

/*
 * -------------------------------------------------
 *  If running in parallel, merge Prodigal results
 * -------------------------------------------------
 */
    if ( params.parallel_split_fasta ) {
            MERGE_PRODIGAL(PRODIGAL.out.amino_acid_fasta.groupTuple(),"faa")

            MERGE_PRODIGAL.out.merged
                .set{merged_prodigal}
    } else {
        PRODIGAL.out.amino_acid_fasta
            .set{merged_prodigal}
    }

/*
 * -------------------------------------------------
 *  OPTIONAL: Run Diamond BLASTp and split data by taxonomy
 * -------------------------------------------------
 */
    if (params.taxonomy_aware) {
        // check if user has given permission for large downloads
        if (params.large_downloads_permission) {
            // Download and prep necessary databases
            PREPARE_NR_DB()
            PREPARE_NR_DB.out.diamond_db
                .set{diamond_db}
            PREPARE_TAXONOMY_DATABASES()
            PREPARE_TAXONOMY_DATABASES.out.taxdump
                .set{ncbi_taxdump}
            PREPARE_TAXONOMY_DATABASES.out.prot_accession2taxid
                .set{prot_accession2taxid}
        } else {
            diamond_db = file("${nr_dmnd_dir}/nr.dmnd")
            ncbi_taxdump = file("${taxdump_tar_gz_dir}/taxdump.tar.gz")
            prot_accession2taxid = file("${prot_accession2taxid_gz_dir}/prot.accession2taxid.gz")
        }

        DIAMOND_BLASTP(merged_prodigal, diamond_db)


        TAXON_ASSIGNMENT(SEQKIT_FILTER.out.fasta, DIAMOND_BLASTP.out.diamond_results, file(taxdump_tar_gz_dir))
        TAXON_ASSIGNMENT.out.taxonomy
            .set{taxonomy_results}

        if (TAXON_ASSIGNMENT.out.bacteria){
            println "ASd"

        } else {
            println "s"
        }
        if (TAXON_ASSIGNMENT.out.archaea){
            println TAXON_ASSIGNMENT.out.archaea

        } else {
            println "s"
        }

    }

    if (params.taxonomy_aware) {
        TAXON_ASSIGNMENT.out.bacteria
        ANALYZE_KMERS(TAXON_ASSIGNMENT.out.bacteria)

    } else {
        ANALYZE_KMERS(SEQKIT_FILTER.out.fasta)
        taxonomy_results = file("$baseDir/assets/dummy_file.txt", checkIfExists: true)
        taxonomy_results = Channel.fromPath( taxonomy_results )

    }

    ANALYZE_KMERS.out.embedded
        .set{kmers_embedded_merged_tsv_ch}

    ANALYZE_KMERS.out.normalized
        .set{kmers_normalized_tsv_ch}

    // --------------------------------------------------------------------------------
    // Run hmmsearch and look for marker genes in contig orfs
    // --------------------------------------------------------------------------------
    MARKERS(PRODIGAL.out.amino_acid_fasta)
 //   HMMER_HMMSEARCH.out.domtblout
 //       .join(PRODIGAL.out.amino_acid_fasta)
 //       .set{hmmsearch_out}
 //   HMMER_HMMSEARCH_FILTER(hmmsearch_out)


    // Before binning we need to merge back everything that was run in parallel
    if ( params.parallel_split_fasta ) {
        MERGE_SPADES_COVERAGE_TSV(SPADES_KMER_COVERAGE.out.coverages.groupTuple(),"coverage")

        MERGE_SPADES_COVERAGE_TSV.out.merged_tsv
            .set{spades_coverage_merged_tsv_ch}

        MERGE_HMMSEARCH(MARKERS.out.markers_tsv.groupTuple(),"markers.tsv")

        MERGE_HMMSEARCH.out.merged_tsv
            .set{markers_tsv_merged_tsv_ch}

    }
    else {
        fasta_ch = SEQKIT_FILTER.out.fasta

        SPADES_KMER_COVERAGE.out.coverages
            .set{spades_coverage_merged_tsv_ch}

        MARKERS.out.markers_tsv
            .set{markers_tsv_merged_tsv_ch}
    }

    BIN_CONTIGS(
        SEQKIT_FILTER.out.fasta,
        kmers_embedded_merged_tsv_ch,
        kmers_normalized_tsv_ch,
        spades_coverage_merged_tsv_ch,
        SEQKIT_FILTER.out.gc_content,
        markers_tsv_merged_tsv_ch,
        taxonomy_results,
        "cluster"
    )

}
