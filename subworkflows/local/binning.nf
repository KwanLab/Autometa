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

        BIN_CONTIGS.out.main
            .join(
                markers
            ).join(
                metagenome
            )
            .set{binning_summary_ch}

        ncbi_tax_dir = file(params.taxdump_tar_gz_dir)

        BINNING_SUMMARY (
            binning_summary_ch,
            binning_column,
            ncbi_tax_dir
        )

    emit:
        binning = BIN_CONTIGS.out.binning
        binning_main = BIN_CONTIGS.out.main
        summary_stats = BINNING_SUMMARY.out.stats
        summary_taxa = BINNING_SUMMARY.out.taxonomies
        metabins = BINNING_SUMMARY.out.metabins
}
