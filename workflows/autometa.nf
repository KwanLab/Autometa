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

include { COVERAGE                    } from '../subworkflows/local/coverage'
include { KMERS                       } from '../subworkflows/local/kmers'
include { PROCESS_METAGENOME          } from '../subworkflows/local/process_metagenome'
include { TAXONOMY_WORKFLOW           } from '../subworkflows/local/taxonomy_workflow'
include { BIN                         } from '../subworkflows/local/binning'

workflow AUTOMETA {

    ch_versions = Channel.empty()

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

    /*
    * -------------------------------------------------
    *  OPTIONAL: Run Diamond BLASTp and split contigs into taxonomic groups
    * -------------------------------------------------
    */

    if (params.taxonomy_aware) {
        TAXONOMY_WORKFLOW (
            filtered_metagenome_fasta,
            orfs_ch
        )
        ch_versions = ch_versions.mix(TAXONOMY_WORKFLOW.out.versions)

        taxonomy_results = TAXONOMY_WORKFLOW.out.taxonomy
        taxdump_files = TAXONOMY_WORKFLOW.out.taxdump_files
        taxonomically_split_fna_ch = TAXONOMY_WORKFLOW.out.taxonomically_split_fna

    } else {
        filtered_metagenome_fasta
            .map { meta, fna ->
                def new_meta = meta.clone()
                new_meta['taxon'] = 'unclassified'
                return [new_meta, fna]
            }
            .set{taxonomically_split_fna_ch}

        Channel
            .fromPath(file("$baseDir/assets/dummy_file.txt", checkIfExists: true ))
            .set{taxonomy_results}
        Channel
            .fromPath(file("$baseDir/assets/dummy_file.txt", checkIfExists: true ))
            .set{taxdump_files}
    }

    /*
    * -------------------------------------------------
    * Calculate k-mer frequencies
    * -------------------------------------------------
    */

    KMERS( taxonomically_split_fna_ch )
    ch_versions = ch_versions.mix(KMERS.out.versions)

    // --------------------------------------------------------------------------------
    // Run hmmscan and look for marker genes in contig orfs
    // --------------------------------------------------------------------------------
    Channel
        .fromList(['bacteria', 'archaea'])
        .set { kingdoms }

    // Ensure orfs_ch is defined before using
    orfs_ch
        .combine(kingdoms)
        .map { pair ->
            def (meta, orfs_file, kingdom) = pair // Correctly extract values from pair
            def new_meta = meta.clone()
            new_meta['taxon'] = kingdom
            return [new_meta, orfs_file]
        }
        .set { orfs_taxon_ch }

    MARKERS( orfs_taxon_ch )

    ch_versions = ch_versions.mix(MARKERS.out.versions)

    markers_ch = MARKERS.out.markers_tsv

    BIN(
        taxonomically_split_fna_ch,
        PROCESS_METAGENOME.out.filtered_metagenome_gc_content,
        markers_ch,
        coverage_ch,
        taxonomy_results,
        KMERS.out.embedded,
        taxdump_files,
        TAXONOMY_WORKFLOW.out.dbtype
    )

    // if (params.mock_test){
    //     BIN.out.binning_results
    //         .join(PROCESS_METAGENOME.out.assembly_to_locus)
    //         .join(PROCESS_METAGENOME.out.assembly_report)
    //         .set { mock_input_ch }

    //     MOCK_DATA_REPORT(
    //         mock_input_ch,
    //         file("$baseDir/lib/mock_data_report.Rmd")
    //     )
    //     ch_versions = ch_versions.mix(MOCK_DATA_REPORT.out.versions)
    // }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}
