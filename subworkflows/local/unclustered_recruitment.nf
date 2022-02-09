params.binning_options                  = [:]
params.unclustered_recruitment_options  = [:]
params.binning_summary_options          = [:]
params.taxdump_tar_gz_dir               = [:]

include { RECRUIT                                         } from './../../modules/local/unclustered_recruitment.nf' addParams( options: params.unclustered_recruitment_options  )
include { BINNING_SUMMARY as UNCLUSTERED_BINNING_SUMMARY  } from './../../modules/local/binning_summary.nf'         addParams( options: params.binning_summary_options, taxdump_tar_gz_dir: params.taxdump_tar_gz_dir )


workflow UNCLUSTERED_RECRUITMENT {

    take:
        metagenome
        kmers_normalized
        coverage
        markers
        taxon_assignments
        binning

    main:

        kmers_normalized
            .join(
                coverage
            ).join(
                binning //BINNING.out.binning
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

        RECRUIT (
            unclustered_recruitment_ch
        )

        RECRUIT.out.main
            .join(
                markers
            ).join(
                metagenome
            )
            .set{unclustered_recruitment_summary_ch}

       // UNCLUSTERED_BINNING_SUMMARY (
       //     unclustered_recruitment_summary_ch,
       //     "recruited_cluster"
       // )

    emit:
        recruitment = RECRUIT.out.binning
        recruitment_main = RECRUIT.out.main
        all_binning_results = binning | mix(RECRUIT.out) | collect
      //  unclustered_recruitment_summary_stats = UNCLUSTERED_BINNING_SUMMARY.out.stats
      //  unclustered_recruitment_summary_taxa = UNCLUSTERED_BINNING_SUMMARY.out.taxonomies
      //  unclustered_recruitment_metabins = UNCLUSTERED_BINNING_SUMMARY.out.metabins
}
