/*
 * -------------------------------------------------
 * Autometa workflow
 * -------------------------------------------------
*/

def modules = params.modules.clone()

def check_for_file(path) {
    return
}

// check if user wants to separate contigs based on taxonomy before binning

if (params.single_db_dir) {
    internal_nr_dmnd_dir = params.single_db_dir
    internal_prot_accession2taxid_gz_dir = params.single_db_dir
    internal_taxdump_tar_gz_dir = params.single_db_dir
}
// TODO: when implementing the ability to set individual DB dirs
// just override e.g. 'internal_nr_dmnd_location' here so users can set
// 'single_db_dir' but also set individual other db paths if they have them
// e.g. if they have nr.dmnd but not the other files.

if (params.large_downloads_permission) {
    // TODO: check if files already exist, if they don't fail the pipeline early at this stage
} else {
    // TODO: check if files exist, if they don't fail the pipeline early at this stage
}


// if these are still null then it meand they weren't set, so make them null.
// this only works because the markov models are inside the docker image.
// that needs to be changed in future versions

if (!params.taxonomy_aware) {
    single_db_dir = null
    internal_nr_dmnd_dir = null
    internal_prot_accession2taxid_gz_dir = null
    internal_taxdump_tar_gz_dir = null
}





/*
 * -------------------------------------------------
 *  Import local modules
 * -------------------------------------------------
*/

include { ANALYZE_KMERS                                         } from '../modules/local/analyze_kmers'               addParams( options: modules['analyze_kmers']                                       )
include { GET_SOFTWARE_VERSIONS                                 } from '../modules/local/get_software_versions'       addParams( options: [publish_files : ['csv':'']]                                   )
include { HMMER_HMMSEARCH                                       } from '../modules/local/hmmer_hmmsearch'            addParams( options: modules['hmmsearch_options']                                   )
include { HMMER_HMMSEARCH_FILTER                                } from '../modules/local/hmmer_hmmsearch_filter'      addParams( options: modules['hmmsearch_filter_options']                            )
include { SEQKIT_FILTER                                         } from '../modules/local/seqkit_filter'               addParams( options: [publish_files : ['*':'']]                                     )
include { MERGE_TSV_WITH_HEADERS as MERGE_SPADES_COVERAGE_TSV   } from '../modules/local/merge_tsv'                   addParams( options: modules['spades_kmer_coverage']                                )
include { MERGE_TSV_WITH_HEADERS as MERGE_HMMSEARCH             } from '../modules/local/merge_tsv'                   addParams( options: modules['merge_hmmsearch_options']                             )
include { SEQKIT_SPLIT                                          } from '../modules/local/seqkit_split'                addParams( options: modules['seqkit_split_options'], num_splits: params.num_splits )
include { SPADES_KMER_COVERAGE                                  } from '../modules/local/spades_kmer_coverage'        addParams( options: modules['spades_kmer_coverage']                                )
include { MERGE_FASTA as MERGE_PRODIGAL                         } from '../modules/local/merge_fasta'                 addParams( )
include { MARKERS                                               } from '../modules/local/markers'                     addParams( options: modules['seqkit_split_options']                                )

/*
 * -------------------------------------------------
 *  Import nf-core modules
 * -------------------------------------------------
*/
// https://github.com/nf-core/modules/tree/master/modules
// https://nf-co.re/tools/#modules
// nf-core modules --help
include { PRODIGAL } from './../modules/nf-core/modules/prodigal/main'  addParams( options: modules['prodigal_options']                         )
/*
 * -------------------------------------------------
 *  Import local subworkflows
 * -------------------------------------------------
*/
include { BIN_CONTIGS      } from '../subworkflows/local/bin_contigs'      addParams( binning_options: modules['binning_options'], unclustered_recruitment_options: modules['unclustered_recruitment_options'], binning_summary_options: modules['binning_summary_options'], prot_accession2taxid_gz_dir: internal_prot_accession2taxid_gz_dir )
include { INPUT_CONTIGS    } from '../subworkflows/local/input_check'      addParams( )
include { CREATE_MOCK      } from '../subworkflows/local/mock_data'        addParams( )
include { TAXON_ASSIGNMENT } from '../subworkflows/local/taxon_assignment' addParams( options: modules['taxon_assignment'], majority_vote_options: modules['majority_vote_options'], split_kingdoms_options: modules['split_kingdoms_options'], nr_dmnd_dir: internal_nr_dmnd_dir, taxdump_tar_gz_dir: internal_taxdump_tar_gz_dir, prot_accession2taxid_gz_dir: internal_prot_accession2taxid_gz_dir, diamond_blastp_options: modules['diamond_blastp_options'], large_downloads_permission: params.large_downloads_permission    )

workflow AUTOMETA {
    ch_software_versions = Channel.empty()

    if (params.mock_test){
        CREATE_MOCK ()
        CREATE_MOCK.out.fasta
            .set{input_ch}
    } else {
        INPUT_CONTIGS ()
        INPUT_CONTIGS.out.metagenome
            .set{input_ch}
    }


    SEQKIT_FILTER(
        input_ch
        )

    // Split contigs FASTA if running in parallel
    if ( params.num_splits > 1 ) {
        SEQKIT_SPLIT (
            SEQKIT_FILTER.out.fasta
        )
        fasta_ch = SEQKIT_SPLIT.out.fasta.transpose()
    } else {
        fasta_ch = SEQKIT_FILTER.out.fasta
    }

/*
 * -------------------------------------------------
 *  Find coverage, currently only pulling from SPADES output
 * -------------------------------------------------
*/

    SPADES_KMER_COVERAGE (
        fasta_ch
    )

/*
 * -------------------------------------------------
 *  Find open reading frames with Prodigal
 * -------------------------------------------------
*/

    PRODIGAL (
        fasta_ch,
        "gbk"
    )

/*
 * -------------------------------------------------
 *  If running in parallel, merge Prodigal results
 * -------------------------------------------------
*/

    if ( params.num_splits > 0 ) {
        MERGE_PRODIGAL (
            PRODIGAL.out.amino_acid_fasta.groupTuple(),
            "faa"
        )
        MERGE_PRODIGAL.out.merged
            .set{merged_prodigal}
    } else {
        PRODIGAL.out.amino_acid_fasta
            .set{merged_prodigal}
    }

/*
 * -------------------------------------------------
 *  OPTIONAL: Run Diamond BLASTp and split contigs into taxonomic groups
 * -------------------------------------------------
*/

    if (params.taxonomy_aware) {
        TAXON_ASSIGNMENT (
            SEQKIT_FILTER.out.fasta,
            merged_prodigal
        )
        TAXON_ASSIGNMENT.out.taxonomy
            .set{taxonomy_results}

        TAXON_ASSIGNMENT.out.bacteria
        ANALYZE_KMERS (
            TAXON_ASSIGNMENT.out.bacteria
            )
    } else {
        ANALYZE_KMERS ( SEQKIT_FILTER.out.fasta )
        taxonomy_results = file( "$baseDir/assets/dummy_file.txt", checkIfExists: true )
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
 // To move to hmmsearch instead of hmmscan:
 // HMMER_HMMSEARCH.out.domtblout
 //       .join(PRODIGAL.out.amino_acid_fasta)
 //       .set{hmmsearch_out}
 // HMMER_HMMSEARCH_FILTER(hmmsearch_out)

    // Before binning we need to merge back everything that was run in parallel
    if ( params.num_splits > 0 ) {
        MERGE_SPADES_COVERAGE_TSV (
            SPADES_KMER_COVERAGE.out.coverages.groupTuple(),
            "coverage"
        )
        MERGE_SPADES_COVERAGE_TSV.out.merged_tsv
            .set{spades_coverage_merged_tsv_ch}

        MERGE_HMMSEARCH (
            MARKERS.out.markers_tsv.groupTuple(),
            "markers.tsv"
        )
        MERGE_HMMSEARCH.out.merged_tsv
            .set{markers_tsv_merged_tsv_ch}
    } else {
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
