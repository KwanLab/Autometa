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
include { PRODIGAL                      } from '../modules/local/prodigal/main.nf'


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

    /*
    * -------------------------------------------------
    *  Init marker sets used for contamination/completion
    *  evalution during binning process
    *  Currently only "bacteria", "archaea"
    *  These are used when taxonomic splitting isn't performed
    * -------------------------------------------------
    */
    taxa_with_marker_sets    = ["bacteria", "archaea"]
    ch_taxa_with_marker_sets = Channel.fromList(taxa_with_marker_sets)

    /*
    * -------------------------------------------------
    *  Process sample list or mock input; filter reads
    * -------------------------------------------------
    */
    PROCESS_METAGENOME()

    filtered_metagenome_fasta = PROCESS_METAGENOME.out.filtered_metagenome_fasta
    ch_versions               = ch_versions.mix(PROCESS_METAGENOME.out.versions)

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
    *  Find open reading frames (proteins) with Prodigal
    * -------------------------------------------------
    */
    PRODIGAL (
        filtered_metagenome_fasta,
        "gbk"
    )

    orfs_ch     = PRODIGAL.out.amino_acid_fasta
    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

    /*
    * -------------------------------------------------
    *  If requested, split contigs into taxonomic groups
    * -------------------------------------------------
    */
    if (params.taxonomy_aware) {

        TAXON_ASSIGNMENT (
            filtered_metagenome_fasta,
            orfs_ch
        )

        ch_taxonomy_results  = TAXON_ASSIGNMENT.out.contig_taxonomy_tsv // [[meta.id, meta.taxon], tsv]
        ch_kmers_input_fasta = TAXON_ASSIGNMENT.out.taxon_split_fasta   // [[meta.id, meta.taxon], fasta]
        ch_taxdump_files     = TAXON_ASSIGNMENT.out.taxdump_files
        ch_versions          = ch_versions.mix(TAXON_ASSIGNMENT.out.versions)

    } else {

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
    if (params.taxonomy_aware) {

        KMERS(ch_kmers_input_fasta )
        ch_kmers_embedded   = KMERS.out.embedded
        ch_kmers_normalized = KMERS.out.normalized
        ch_versions         = ch_versions.mix(KMERS.out.versions)





    } else {
        // prevent computing kmers multiple times on the same data
        KMERS(filtered_metagenome_fasta )
        ch_versions         = ch_versions.mix(KMERS.out.versions)

        KMERS.out.embedded
            .combine(
                Channel
                    .fromList(
                        taxa_with_marker_sets
                        )
            )
            .map{ meta, fasta, taxon ->
                    [meta + [taxon: taxon], fasta]
                }
                .set{ ch_kmers_embedded }

        KMERS.out.normalized
            .combine(
                Channel
                    .fromList(
                        taxa_with_marker_sets
                        )
            )
            .map{ meta, fasta, taxon ->
                    [meta + [taxon: taxon], fasta]
                }
                .set{ ch_kmers_normalized }

        // This adds taxon to the TAXON_ASSIGNMENT.out.taxonomy meta map
        // e.g. [[sample_id, taxon:bacteria], fastapath]
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
    }

    /*
    * -------------------------------------------------
    * Map channels and taxonomic
    * -------------------------------------------------
    */
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

    orfs_ch
        .combine(ch_taxa_with_marker_sets)
        .map{ meta, fasta, taxon ->
                [meta + [taxon: taxon], fasta]
            }
        .set{ ch_new_orf }

    /*
    * -------------------------------------------------
    * Run hmmscan and look for marker genes in contig orfs
    * -------------------------------------------------
    */
    MARKERS(ch_new_orf)
    ch_markers  = MARKERS.out.markers_tsv
    ch_versions = ch_versions.mix(MARKERS.out.versions)

    /*
    * -------------------------------------------------
    * Prepare inputs for binning channel
    * -------------------------------------------------
    */
    ch_kmers_embedded
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

    /*
    * -------------------------------------------------
    * Bin contigs
    * -------------------------------------------------
    */
    BINNING(
        ch_binning
    )
    ch_versions = ch_versions.mix(BINNING.out.versions)

    /*
    * -------------------------------------------------
    * Recruit unclustered contigs
    * -------------------------------------------------
    */

    if (params.unclustered_recruitment) {
        // Prepare inputs for recruitment channel
        ch_kmers_normalized
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

    /*
    * -------------------------------------------------
    * Summarize binning results
    * -------------------------------------------------
    */
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

    /*
    * -------------------------------------------------
    * Create reports for test runs (testing the workflow)
    * -------------------------------------------------
    */

    if (workflow.profile.contains("test")){

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

    /*
    * -------------------------------------------------
    * Collect software versions from all processes that
    * ran into a single output yeahml! file
    * -------------------------------------------------
    */
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}
