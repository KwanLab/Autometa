/*
 * -------------------------------------------------
 * Autometa workflow
 * -------------------------------------------------
*/

def modules = params.modules.clone()

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

// if these are still null then it means they weren't set, so make them null.
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

include { ANALYZE_KMERS                                         } from '../modules/local/analyze_kmers'               addParams( options: modules['analyze_kmers_options'] )
include { GET_SOFTWARE_VERSIONS                                 } from '../modules/local/get_software_versions'       addParams( options: [publish_files : ['csv':'']]     )
include { SEQKIT_FILTER                                         } from '../modules/local/seqkit_filter'               addParams( options: [publish_files : ['*':'']]       )
include { SPADES_KMER_COVERAGE                                  } from '../modules/local/spades_kmer_coverage'        addParams( options: modules['spades_kmer_coverage']  )
include { MARKERS                                               } from '../modules/local/markers'                     addParams( options: modules['seqkit_split_options']  )
include { BIN_CONTIGS                                           } from '../modules/local/bin_contigs'                 addParams( options: modules['bin_contigs_options']   )
include { RECRUIT                                               } from '../modules/local/unclustered_recruitment'     addParams( options: modules['unclustered_recruitment_options'])
include { BINNING_SUMMARY                                       } from '../modules/local/binning_summary'             addParams( options: modules['binning_summary_options']   )
include { MOCK_DATA_REPORT                                      } from '../modules/local/mock_data_reporter'          addParams( options: modules['mock_data_report']      )

/*
 * -------------------------------------------------
 *  Import nf-core modules
 * -------------------------------------------------
*/
// https://github.com/nf-core/modules/tree/master/modules
// https://nf-co.re/tools/#modules
// nf-core modules --help
include { PRODIGAL } from './../modules/nf-core/modules/prodigal/main'  addParams( options: modules['prodigal_options'] )

/*
 * -------------------------------------------------
 *  Import local subworkflows
 * -------------------------------------------------
*/

include { CREATE_MOCK              } from '../subworkflows/local/mock_data'                addParams( get_genomes_for_mock: modules['get_genomes_for_mock'])
include { INPUT_CHECK              } from '../subworkflows/local/input_check'              addParams( )
include { CONTIG_COVERAGE          } from '../subworkflows/local/contig_coverage'          addParams( align_reads_options: modules['align_reads_options'], samtools_viewsort_options: modules['samtools_viewsort_options'], bedtools_genomecov_options: modules['bedtools_genomecov_options'])
include { TAXON_ASSIGNMENT         } from '../subworkflows/local/taxon_assignment'         addParams( options: modules['taxon_assignment'], majority_vote_options: modules['majority_vote_options'], split_kingdoms_options: modules['split_kingdoms_options'], nr_dmnd_dir: internal_nr_dmnd_dir, taxdump_tar_gz_dir: internal_taxdump_tar_gz_dir, prot_accession2taxid_gz_dir: internal_prot_accession2taxid_gz_dir, diamond_blastp_options: modules['diamond_blastp_options'], large_downloads_permission: params.large_downloads_permission )

workflow AUTOMETA {
    ch_software_versions = Channel.empty()
    samplesheet_ch = Channel.fromPath(params.input)

    if (params.mock_test){
        CREATE_MOCK()
        CREATE_MOCK.out.fasta
            .set{metagenome_ch}
        coverage_tab_ch = Channel.empty()
    } else {
        INPUT_CHECK(samplesheet_ch)
        INPUT_CHECK.out.metagenome
            .set{metagenome_ch}
        INPUT_CHECK.out.coverage
            .set{coverage_tab_ch}
    }


    SEQKIT_FILTER(
        metagenome_ch
    )
    fasta_ch = SEQKIT_FILTER.out.fasta

    /*
    * -------------------------------------------------
    *  Find coverage, currently only pulling from SPADES output
    * -------------------------------------------------
    */


    if (!params.mock_test) {
        coverage_input_ch = fasta_ch.join(INPUT_CHECK.out.reads)
    } else {
        coverage_input_ch = Channel.empty()
    }

    CONTIG_COVERAGE (
        coverage_input_ch
    )
    CONTIG_COVERAGE.out.coverage
        .set{contig_coverage_ch}

    SPADES_KMER_COVERAGE (
        fasta_ch,
    )
    SPADES_KMER_COVERAGE.out.coverage
        .set{spades_kmer_coverage_ch}
    // https://nextflow-io.github.io/patterns/index.html#_conditional_process_executions
    contig_coverage_ch
        .mix(spades_kmer_coverage_ch)
        .mix(coverage_tab_ch)
        .set{coverage_ch}

    /*
    * -------------------------------------------------
    *  Find open reading frames with Prodigal
    * -------------------------------------------------
    */

    PRODIGAL (
        fasta_ch,
        "gbk"
    )

    PRODIGAL.out.amino_acid_fasta
        .set{orfs_ch}

    /*
    * -------------------------------------------------
    *  OPTIONAL: Run Diamond BLASTp and split contigs into taxonomic groups
    * -------------------------------------------------
    */

    if (params.taxonomy_aware) {
        TAXON_ASSIGNMENT (
            fasta_ch,
            orfs_ch
        )
        TAXON_ASSIGNMENT.out.taxonomy
            .set{taxonomy_results}
        if (params.kingdom.equals('bacteria')) {
            kmers_input_ch = TAXON_ASSIGNMENT.out.bacteria
        } else {
            // params.kingdom.equals('archaea')
            kmers_input_ch = TAXON_ASSIGNMENT.out.archaea
        }

    } else {
        kmers_input_ch = fasta_ch
        taxonomy_results = file( "$baseDir/assets/dummy_file.txt", checkIfExists: true )
        taxonomy_results = Channel.fromPath( taxonomy_results )
    }

    /*
    * -------------------------------------------------
    * Calculate k-mer frequencies
    * -------------------------------------------------
    */

    ANALYZE_KMERS(
        kmers_input_ch
    )
    ANALYZE_KMERS.out.normalized
        .set{kmers_normalized_ch}

    ANALYZE_KMERS.out.embedded
        .set{kmers_embedded_ch}


    // --------------------------------------------------------------------------------
    // Run hmmscan and look for marker genes in contig orfs
    // --------------------------------------------------------------------------------

    MARKERS(
        orfs_ch
    )
    MARKERS.out.markers_tsv
        .set{markers_ch}

    // Prepare inputs for binning channel
    kmers_embedded_ch
        .join(coverage_ch)
        .join(SEQKIT_FILTER.out.gc_content)
        .join(markers_ch)
        .set{binning_ch}
    if (params.taxonomy_aware) {
        binning_ch
            .join(taxonomy_results)
            .set{binning_ch}
    } else {
        binning_ch
            .combine(taxonomy_results)
            .set{binning_ch}
    }

    BIN_CONTIGS(
        binning_ch
    )

    if (params.unclustered_recruitment) {
        // Prepare inputs for recruitment channel
        kmers_normalized_ch
            .join(coverage_ch)
            .join(BIN_CONTIGS.out.binning)
            .join(markers_ch)
            .set{recruitment_ch}
        if (params.taxonomy_aware) {
            recruitment_ch
                .join(taxonomy_results)
                .set{recruitment_ch}
        } else {
            recruitment_ch
                .combine(taxonomy_results)
                .set{recruitment_ch}
        }
        RECRUIT(
            recruitment_ch
        )
        RECRUIT.out.main
            .set{binning_results_ch}
        binning_col = Channel.value("recruited_cluster")
    } else {
        BIN_CONTIGS.out.main
            .set{binning_results_ch}
        binning_col = Channel.value("cluster")
    }

    // Set inputs for binning summary
    binning_results_ch
        .join(markers_ch)
        .join(fasta_ch)
        .set{binning_summary_ch}

    ncbi = file(params.single_db_dir)

    BINNING_SUMMARY(
        binning_summary_ch,
        binning_col,
        ncbi,
    )

    if (params.mock_test){
        binning_results_ch
            .join(CREATE_MOCK.out.assembly_to_locus)
            .join(CREATE_MOCK.out.assembly_report)
            .set { mock_input_ch }

        MOCK_DATA_REPORT(
            mock_input_ch,
            file("$baseDir/lib/mock_data_report.Rmd")
        )
    }

}
