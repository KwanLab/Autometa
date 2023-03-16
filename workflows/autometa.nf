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
include { CUSTOM_DUMPSOFTWAREVERSIONS             } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { MARKERS                                 } from '../modules/local/markers'
include { BINNING                                 } from '../modules/local/binning'
include { RECRUIT                                 } from '../modules/local/unclustered_recruitment'
include { BINNING_SUMMARY                         } from '../modules/local/binning_summary'
include { MOCK_DATA_REPORT                        } from '../modules/local/mock_data_reporter'

/*
 * -------------------------------------------------
 *  Import nf-core modules
 * -------------------------------------------------
*/
// https://github.com/nf-core/modules/tree/master/modules
// https://nf-co.re/tools/#modules
// nf-core modules --help
include { PRODIGAL } from './../modules/nf-core/prodigal/main.nf'

/*
 * -------------------------------------------------
 *  Import local subworkflows
 * -------------------------------------------------
*/

include { COVERAGE                       } from '../subworkflows/local/coverage'
include { KMERS                       } from '../subworkflows/local/kmers'
include { PROCESS_METAGENOME          } from '../subworkflows/local/process_metagenome'
include { TAXON_ASSIGNMENT            } from '../subworkflows/local/taxon_assignment'

workflow AUTOMETA {

    ch_versions = Channel.empty()

    // only these taxa can be binned
    taxa_with_marker_sets = ["bacteria", "archaea"]
    taxa_with_marker_sets_ch = Channel.fromList(taxa_with_marker_sets)

    PROCESS_METAGENOME()
    ch_versions = ch_versions.mix(PROCESS_METAGENOME.out.versions)

    COVERAGE(
        PROCESS_METAGENOME.out.filtered_metagenome_fasta,
        PROCESS_METAGENOME.out.filtered_metagenome_fasta_and_reads,
        PROCESS_METAGENOME.out.user_provided_coverage_table
    )
    ch_versions = ch_versions.mix(COVERAGE.out.versions)

    filtered_metagenome_fasta = PROCESS_METAGENOME.out.filtered_metagenome_fasta
    coverage_ch = COVERAGE.out.coverage_ch

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

        taxonomy_results_ch     = TAXON_ASSIGNMENT.out.contig_taxonomy_tsv 
        kmers_input_fasta_ch    = TAXON_ASSIGNMENT.out.taxon_split_fasta
        found_taxa_list_ch      = TAXON_ASSIGNMENT.out.found_taxa_list
        taxdump_files_ch        = TAXON_ASSIGNMENT.out.taxdump_files

    } else {

        found_taxa_list_ch = taxa_with_marker_sets_ch
            
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
                .set{ kmers_input_fasta_ch }

        Channel
            .fromPath(file("$baseDir/assets/dummy_file.txt", checkIfExists: true ))
            .set{taxonomy_results_ch}
        Channel
            .fromPath(file("$baseDir/assets/dummy_file.txt", checkIfExists: true ))
            .set{taxdump_files_ch}

    }

    /*
    * -------------------------------------------------
    * Calculate k-mer frequencies
    * -------------------------------------------------
    */

    KMERS(kmers_input_fasta_ch )
    ch_versions = ch_versions.mix(KMERS.out.versions)

    // --------------------------------------------------------------------------------
    // Run hmmscan and look for marker genes in contig orfs
    // --------------------------------------------------------------------------------

    orfs_ch
        .combine(taxa_with_marker_sets_ch)
        .map{ meta, fasta, taxon ->                    
                [meta + [taxon: taxon], fasta]
            }
        .set{ new_orf_ch }

    MARKERS(new_orf_ch)
    ch_versions = ch_versions.mix(MARKERS.out.versions)

    markers_ch = MARKERS.out.markers_tsv

    coverage_ch
        .combine(taxa_with_marker_sets_ch)
        .map{ meta, fasta, taxon ->                    
                [meta + [taxon: taxon], fasta]
            }
        .set{ new_coverage_ch }

    PROCESS_METAGENOME.out.filtered_metagenome_gc_content
        .combine(taxa_with_marker_sets_ch)
        .map{ meta, fasta, taxon ->                    
                [meta + [taxon: taxon], fasta]
            }
        .set{ new_filtered_metagenome_gc_content_ch }

    // markers_ch.map{k,v -> println "$k and $v"}

    // Prepare inputs for binning channel
    KMERS.out.embedded
        .join(new_coverage_ch)
        .join(new_filtered_metagenome_gc_content_ch)
        .join(markers_ch)
        .set{binning_ch}



    if (params.taxonomy_aware) {
        binning_ch
            .join(taxonomy_results_ch)
            .set{binning_ch}
    } else {
        binning_ch
            .combine(taxonomy_results_ch)
            .set{binning_ch}
    }
    



    BINNING(
        binning_ch
    )

}
