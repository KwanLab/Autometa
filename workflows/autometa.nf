/*
 * -------------------------------------------------
 * Autometa workflow
 * -------------------------------------------------
*/

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
        TAXON_ASSIGNMENT (
            filtered_metagenome_fasta,
            orfs_ch
        )
        ch_versions = ch_versions.mix(TAXON_ASSIGNMENT.out.versions)

        taxonomy_results = TAXON_ASSIGNMENT.out.taxonomy
        taxdump_files = TAXON_ASSIGNMENT.out.taxdump_files

        if (params.kingdom.equals('bacteria')) {
            kmers_input_ch = TAXON_ASSIGNMENT.out.bacteria
        } else {
            // params.kingdom.equals('archaea')
            kmers_input_ch = TAXON_ASSIGNMENT.out.archaea
        }
    } else {
        kmers_input_ch = filtered_metagenome_fasta
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

    KMERS( kmers_input_ch )
    ch_versions = ch_versions.mix(KMERS.out.versions)

    // --------------------------------------------------------------------------------
    // Run hmmscan and look for marker genes in contig orfs
    // --------------------------------------------------------------------------------

    MARKERS( orfs_ch )
    ch_versions = ch_versions.mix(MARKERS.out.versions)

    markers_ch = MARKERS.out.markers_tsv

    // Prepare inputs for binning channel
    KMERS.out.embedded
        .join(coverage_ch)
        .join(PROCESS_METAGENOME.out.filtered_metagenome_gc_content)
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

    BINNING(
        binning_ch
    )
    ch_versions = ch_versions.mix(BINNING.out.versions)

    if (params.unclustered_recruitment) {
        // Prepare inputs for recruitment channel
        KMERS.out.normalized
            .join(coverage_ch)
            .join(BINNING.out.main)
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
        ch_versions = ch_versions.mix(RECRUIT.out.versions)

        RECRUIT.out.main
            .set{binning_results_ch}
        binning_col = Channel.of("recruited_cluster")
    } else {
        binning_results_ch = BINNING.out.main
        binning_col = Channel.of("cluster")
    }

    // Set inputs for binning summary
    binning_results_ch
        .join(markers_ch)
        .join(filtered_metagenome_fasta)
        .combine(binning_col)
        .set{binning_summary_ch}


    BINNING_SUMMARY(
        binning_summary_ch,
        taxdump_files.collect()
    )



    ch_versions = ch_versions.mix(BINNING_SUMMARY.out.versions)

    if (params.mock_test){
        binning_results_ch
            .join(PROCESS_METAGENOME.out.assembly_to_locus)
            .join(PROCESS_METAGENOME.out.assembly_report)
            .set { mock_input_ch }

        MOCK_DATA_REPORT(
            mock_input_ch,
            file("$baseDir/lib/mock_data_report.Rmd")
        )
        ch_versions = ch_versions.mix(MOCK_DATA_REPORT.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

}
