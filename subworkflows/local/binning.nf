params.binning_options                  = [:]
params.binning_summary_options          = [:]
params.taxdump_tar_gz_dir               = [:]

include { BIN_CONTIGS     } from './../../modules/local/bin_contigs.nf'             addParams( options: params.binning_options                                                        )
include { BINNING_SUMMARY } from './../../modules/local/binning_summary.nf'         addParams( options: params.binning_summary_options, taxdump_tar_gz_dir: params.taxdump_tar_gz_dir )


workflow BINNING {

    take:
        metagenome
        kmers_embedded
        coverage
        gc_content
        markers
        taxon_assignments
        binning_column

    main:
        kmers_embedded
            .join(
                coverage
                )
            .join(
                gc_content
                )
            .join(
                markers
                )
            .set{metagenome_annotations}

        if (params.taxonomy_aware) {
            metagenome_annotations
                .join(
                    taxon_assignments
                )
                .set{binning_ch}
        } else {
            metagenome_annotations
                .combine(
                    taxon_assignments
                )
                .set{binning_ch}
        }

        BIN_CONTIGS (
            binning_ch
        )

        kmers_embedded
            .join(
                coverage
            ).join(
                BIN_CONTIGS.out.binning
            ).join(
                markers
            )
            .set{coverage_binningout_markers}

        if (params.taxonomy_aware) {
            coverage_binningout_markers
                .join(
                    taxon_assignments
                )
                .set{unclustered_recruitment_ch}
        } else {
            coverage_binningout_markers
                .combine(
                    taxon_assignments
                )
                .set{unclustered_recruitment_ch}
        }

        BIN_CONTIGS.out.main
            .join(
                markers
            ).join(
                metagenome
            )
            .set{binning_summary_ch}

        BINNING_SUMMARY (
            binning_summary_ch,
            binning_column
        )

    emit:
        binning = BIN_CONTIGS.out.binning
        binning_main = BIN_CONTIGS.out.main
        summary_stats = BINNING_SUMMARY.out.stats
        summary_taxa = BINNING_SUMMARY.out.taxonomies
        metabins = BINNING_SUMMARY.out.metabins
}
