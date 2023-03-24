/*
 * -------------------------------------------------
 * Autometa workflow
 * -------------------------------------------------
*/

/*
 * -------------------------------------------------
 *  Import local modules
 * -------------------------------------------------
*/
include { CUSTOM_DUMPSOFTWAREVERSIONS   } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MARKERS                       } from '../modules/local/markers'
include { BINNING                       } from '../modules/local/binning'
include { UNCLUSTERED_RECRUIT           } from '../modules/local/unclustered_recruitment'
include { BINNING_SUMMARY               } from '../modules/local/binning_summary'
include { MOCK_DATA_REPORT              } from '../modules/local/mock_data_reporter'

/*
 * -------------------------------------------------
 *  Import nf-core modules
 * -------------------------------------------------
*/
// https://github.com/nf-core/modules/tree/master/modules
// https://nf-co.re/tools/#modules
// nf-core modules --help
include { PRODIGAL                      } from './../modules/nf-core/prodigal/main.nf'

/*
 * -------------------------------------------------
 *  Import local subworkflows
 * -------------------------------------------------
*/

include { COVERAGE                      } from '../subworkflows/local/coverage'
include { KMERS                         } from '../subworkflows/local/kmers'
include { PROCESS_METAGENOME            } from '../subworkflows/local/process_metagenome'
include { TAXON_ASSIGNMENT              } from '../subworkflows/local/taxon_assignment'

workflow AUTOMETA {

    ch_versions = Channel.empty()

    // CUrrently only "bacteria", "archaea" and have markers and can be binned
    taxa_with_marker_sets = ["bacteria", "archaea"]

    ch_taxa_with_marker_sets = Channel.fromList(taxa_with_marker_sets)

    /*
    * -------------------------------------------------
    *  Process sample list or mock input; filter reads
    * -------------------------------------------------
    */
    PROCESS_METAGENOME()

    filtered_metagenome_fasta = PROCESS_METAGENOME.out.filtered_metagenome_fasta
    ch_versions = ch_versions.mix(PROCESS_METAGENOME.out.versions)

    /*
    * -------------------------------------------------
    *  Extract SPADES coverage or align reads to metagenome(s)
    * -------------------------------------------------
    */
    COVERAGE(
        PROCESS_METAGENOME.out.filtered_metagenome_fasta,
        PROCESS_METAGENOME.out.filtered_metagenome_fasta_and_reads,
        PROCESS_METAGENOME.out.user_provided_coverage_table
    )

    ch_coverage = COVERAGE.out.ch_coverage
    ch_versions = ch_versions.mix(COVERAGE.out.versions)

    /*
    * -------------------------------------------------
    *  Find open reading frames with Prodigal
    * -------------------------------------------------
    */

    PRODIGAL (
        filtered_metagenome_fasta,
        "gbk"
    )
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    PRODIGAL.out.amino_acid_fasta
        .set{orfs_ch}

    if (params.taxonomy_aware) {
        /*
        * -------------------------------------------------
        *  OPTIONAL: Run Diamond BLASTp and split contigs into taxonomic groups
        * -------------------------------------------------
        */
        TAXON_ASSIGNMENT (
            filtered_metagenome_fasta,
            orfs_ch
        )
        ch_versions = ch_versions.mix(TAXON_ASSIGNMENT.out.versions)

        ch_taxonomy_results     = TAXON_ASSIGNMENT.out.contig_taxonomy_tsv
        ch_kmers_input_fasta    = TAXON_ASSIGNMENT.out.taxon_split_fasta
        ch_found_taxa_list      = TAXON_ASSIGNMENT.out.found_taxa_list
        ch_taxdump_files        = TAXON_ASSIGNMENT.out.taxdump_files

    } else {

        ch_found_taxa_list = ch_taxa_with_marker_sets

        // This adds taxon to the TAXON_ASSIGNMENT.out.taxonomy meta map
        filtered_metagenome_fasta
            .combine(
                Channel
                    .fromList(
                        taxa_with_marker_sets
                        )
            )
            .map{ meta, fasta, taxon ->
                    [meta + [taxon: taxon], fasta]
                }
                .set{ ch_kmers_input_fasta }

        Channel
            .fromPath(file("$baseDir/assets/dummy_file.txt", checkIfExists: true ))
            .set{ch_taxonomy_results}
        Channel
            .fromPath(file("$baseDir/assets/dummy_file.txt", checkIfExists: true ))
            .set{ch_taxdump_files}

    }

    /*
    * -------------------------------------------------
    * Calculate k-mer frequencies
    * -------------------------------------------------
    */

    KMERS(ch_kmers_input_fasta )
    ch_versions = ch_versions.mix(KMERS.out.versions)

    // --------------------------------------------------------------------------------
    // Run hmmscan and look for marker genes in contig orfs
    // --------------------------------------------------------------------------------

    orfs_ch
        .combine(ch_taxa_with_marker_sets)
        .map{ meta, fasta, taxon ->
                [meta + [taxon: taxon], fasta]
            }
        .set{ ch_new_orf }

    MARKERS(ch_new_orf)
    ch_versions = ch_versions.mix(MARKERS.out.versions)

    ch_markers = MARKERS.out.markers_tsv

    ch_coverage
        .combine(ch_taxa_with_marker_sets)
        .map{ meta, fasta, taxon ->
                [meta + [taxon: taxon], fasta]
            }
        .set{ new_ch_coverage }

    PROCESS_METAGENOME.out.filtered_metagenome_gc_content
        .combine(ch_taxa_with_marker_sets)
        .map{ meta, fasta, taxon ->
                [meta + [taxon: taxon], fasta]
            }
        .set{ ch_new_filtered_metagenome_gc_content }

    // Prepare inputs for binning channel
    KMERS.out.embedded
        .join(new_ch_coverage)
        .join(ch_new_filtered_metagenome_gc_content)
        .join(ch_markers)
        .set{ch_binning}



    if (params.taxonomy_aware) {
        ch_binning
            .join(ch_taxonomy_results)
            .set{ch_binning}
    } else {
        ch_binning
            .combine(ch_taxonomy_results)
            .set{ch_binning}
    }




    BINNING(
        ch_binning
    )
    ch_versions = ch_versions.mix(BINNING.out.versions)

    if (params.unclustered_recruitment) {
        // Prepare inputs for recruitment channel
        KMERS.out.normalized
            .join(new_ch_coverage)
            .join(BINNING.out.main)
            .join(ch_markers)
            .set{recruitment_ch}

        if (params.taxonomy_aware) {
            recruitment_ch
                .join(ch_taxonomy_results)
                .set{recruitment_ch}
        } else {
            recruitment_ch
                .combine(ch_taxonomy_results)
                .set{recruitment_ch}
        }
        UNCLUSTERED_RECRUIT(
            recruitment_ch
        )
        ch_versions = ch_versions.mix(UNCLUSTERED_RECRUIT.out.versions)

        UNCLUSTERED_RECRUIT.out.main
            .set{ch_binning_results}
        binning_col = Channel.of("recruited_cluster")
    } else {
        ch_binning_results = BINNING.out.main
        binning_col = Channel.of("cluster")
    }

    // Set inputs for binning summary
    ch_binning_results
        .join(ch_markers)
        .join(ch_kmers_input_fasta)
        .combine(binning_col)
        .set{binning_summary_ch}

    BINNING_SUMMARY(
        binning_summary_ch,
        ch_taxdump_files
    )
    ch_versions = ch_versions.mix(BINNING_SUMMARY.out.versions)

    if (workflow.stubRun){

        PROCESS_METAGENOME.out.assembly_to_locus
            .combine(ch_taxa_with_marker_sets)
            .map{ meta, fasta, taxon ->
                    [meta + [taxon: taxon], fasta]
                }
            .set{ ch_new_assembly_to_locus }

        PROCESS_METAGENOME.out.assembly_report
            .combine(ch_taxa_with_marker_sets)
            .map{ meta, fasta, taxon ->
                    [meta + [taxon: taxon], fasta]
                }
            .set{ ch_new_assembly_report }

        ch_binning_results
            .join(ch_new_assembly_to_locus)
            .join(ch_new_assembly_report)
            .set { ch_mock_input }

        MOCK_DATA_REPORT(
            ch_mock_input,
            file("$baseDir/lib/mock_data_report.Rmd")
        )
        ch_versions = ch_versions.mix(MOCK_DATA_REPORT.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}
