include { BINNING           } from '../../modules/local/binning'
include { RECRUIT           } from '../../modules/local/unclustered_recruitment'
include { BINNING_SUMMARY   } from '../../modules/local/binning_summary'



workflow BIN {
    take:
        filtered_metagenome_fasta
        filtered_metagenome_gc_content
        markers_ch
        coverage_ch
        taxonomy_results_ch
        taxonomically_split_fna_ch
        taxdump_files
        dbtype

    main:
        ch_versions = Channel.empty()

        // has taxon:
        //     taxonomically_split_fna_ch
        //     markers_ch
        // not has taxon:
        //     coverage_ch
        //     filtered_metagenome_gc_content
        //     taxonomy_results_ch

        // Transform taxonomic-specific channels (keep taxon info)
        taxonomically_split_fna_ch = taxonomically_split_fna_ch
            .map { meta, files ->
                def key = [id: meta.id, taxon: meta.taxon]
                [key, files]
            }

        markers_ch = markers_ch
            .map { meta, files ->
                def key = [id: meta.id, taxon: meta.taxon]
                [key, files]
            }

        // Transform per-sample channels
        coverage_ch = coverage_ch
            .map { meta, files ->
                [meta.id, files]
            }

        filtered_metagenome_gc_content = filtered_metagenome_gc_content
            .map { meta, files ->
                [meta.id, files]
            }

        taxonomy_results_ch = taxonomy_results_ch
            .map { meta, files ->
                [meta.id, files]
            }

        // Create branched workflow
        workflow_branch = taxonomically_split_fna_ch
            .join(markers_ch)
            .map { key, kmers, markers ->
                // Use the full sample ID as the join key while preserving taxon info
                [key.id, [id: key.id, taxon: key.taxon], kmers, markers]
            }
            .combine(coverage_ch, by: 0)
            .combine(filtered_metagenome_gc_content, by: 0)
            .combine(taxonomy_results_ch, by: 0)
            .map { id, meta, kmers, markers, coverage, gc_content, taxonomy_results ->
                // Final structure: [meta with taxon, files...]
                [meta, kmers, markers, coverage, gc_content, taxonomy_results]
            }

        // Set the output channel
        workflow_branch.set { to_bin_ch }

        BINNING(
            to_bin_ch
        )

        ch_versions = ch_versions.mix(BINNING.out.versions)

        if (params.unclustered_recruitment) {
            // Prepare inputs for recruitment channel

            to_bin_ch
                .join(BINNING.out.main)
                .set{recruitment_ch}

            RECRUIT(
                recruitment_ch
            )
            ch_versions = ch_versions.mix(RECRUIT.out.versions)

            RECRUIT.out.main
                .set{binning_results_ch}
            binning_col = Channel.from("recruited_cluster")
        } else {
            binning_results_ch = BINNING.out.main
            binning_col = Channel.from("cluster")
        }


        // Set inputs for binning summary
        binning_results_ch
            .map { meta, files -> [meta.subMap(['id']), meta, files] }
            .join(markers_ch.map { meta, files -> [meta.subMap(['id']), files] })
            .join(filtered_metagenome_fasta.map { meta, files -> [meta.subMap(['id']), files] })
            .map { it.drop(1) }
            .set{binning_summary_input_ch}

        if (params.debug) {
            binning_results_ch.view { meta ->
                println "binning_results_ch: ${meta}"
            }
            markers_ch.view { meta ->
                println "markers_ch: ${meta}"
            }
            filtered_metagenome_fasta.view { meta ->
                println "filtered_metagenome_fasta: ${meta}"
            }
            taxdump_files.view { meta ->
                println "taxdump_files: ${meta}"
            }
            markers_ch.view { meta ->
                println "markers_ch: ${meta}"
            }
            binning_col.view { meta ->
                println "binning_col: ${meta}"
            }
            binning_summary_input_ch.view { meta ->
                println "binning_summary_input_ch: ${meta}"
            }
        }

        binning_summary_input_ch
            .combine(taxdump_files.toList())
            .combine(dbtype)
            .combine(binning_col)
            .set{binning_summary_input_ch2}

        BINNING_SUMMARY(
            binning_summary_input_ch2
        )
        ch_versions = ch_versions.mix(BINNING_SUMMARY.out.versions)

    emit:
        binning_results = binning_results_ch
        // binning_summary = BINNING_SUMMARY
        versions = ch_versions
}
