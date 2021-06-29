
params.binning_options = [:]
params.unclustered_recruitment_options = [:]
params.binning_summary_options      = [:]


include { BINNING } from './../../modules/local/binning.nf' addParams( options: params.binning_options )
include { UNCLUSTERED_RECRUITMENT } from './../../modules/local/unclustered_recruitment.nf' addParams( options: params.unclustered_recruitment_options )
include { BINNING_SUMMARY } from './../../modules/local/binning_summary.nf' addParams( options: params.binning_summary_options )


workflow BIN_CONTIGS {

    take:
    metagenome
    kmers_embedded
    kmers_normalized
    coverage
    gc_content
    markers
    taxon_assignments
    binning_column
    

    main: 
    
        kmers_embedded.join(
            coverage
        ).join(
            gc_content
        ).join(
            markers
        )
        .set{binning_ch}
  
        BINNING(binning_ch, taxon_assignments)
    
        kmers_normalized.join(
            coverage
        ).join(
            BINNING.out.binning
        ).join(
            markers
        )
        .set{unclustered_recruitment_ch}

        UNCLUSTERED_RECRUITMENT(unclustered_recruitment_ch, taxon_assignments)
    
        BINNING.out.main.join(
            markers
        ).join(
            metagenome
        )
        .set{binning_summary_ch}

        BINNING_SUMMARY(binning_summary_ch, binning_column)
  
    emit:
        binning = BINNING.out.binning
        binning_main = BINNING.out.main
        recruitment = UNCLUSTERED_RECRUITMENT.out.binning
        recruitment_main = UNCLUSTERED_RECRUITMENT.out.main
        all_binning_results = BINNING.out.binning | mix(UNCLUSTERED_RECRUITMENT.out) | collect
        summary_stats = BINNING_SUMMARY.out.stats
        summary_taxa = BINNING_SUMMARY.out.taxonomies
        metabins = BINNING_SUMMARY.out.metabins

}
