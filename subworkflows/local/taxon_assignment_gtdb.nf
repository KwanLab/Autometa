
include { PREPARE_GTDB_DB   } from './prepare_gtdb.nf'
include { TAXON_SPLIT       } from './taxon_split.nf'

process EXTRACT_ORFS {
    tag "Extracting ORFs from taxon-assigned metagenome contigs"
    label 'process_low'
    
    conda "bioconda::autometa"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/autometa:2.2.0--pyh7cba7a3_0"
    } else {
        container "jasonkwan/autometa:${params.autometa_image_tag}"
    }

    input:
        tuple val(meta), path(contigs), path (orfs)

    output:
        tuple val(meta), path("${meta.id}.gtdb_input.fna"), path("${meta.id}_gtdb_input_orfs.faa.gz"), emit: split_orfs       
        path "versions.yml", emit: versions
    
    script:
        def prefix = task.ext.prefix ?: "${meta.id}"
        """
        grep -h ">" $contigs | \\
                    sed 's/^>//' | \\
                    cut -f1 -d" " | \\
                    sed 's/\\\$/_/' | \\
                grep -f - $orfs |\\
                    cut -f1 -d" " |\\
                        sed 's/^>//'  > orf_ids
                        
                # Retrieve ORF seqs from ORF IDs
                seqkit grep \
                    --pattern-file orf_ids \
                    --out-file ${meta.id}_gtdb_input_orfs.faa.gz \
                    $orfs

        cat $contigs > ${meta.id}.gtdb_input.fna

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            autometa: \$(autometa --version | sed -e 's/autometa: //g')
            gtdb: $params.gtdb_version
        END_VERSIONS
        """
}


// Autometa taxon assignment workflow
workflow GTDB_TAXON_ASSIGNMENT {
    take:
        split_metagenome_contigs
        orfs

    main:
        ch_versions = Channel.empty()
        dbtype_ch = channel.value( 'gtdb')

        PREPARE_GTDB_DB()
        ch_versions = ch_versions.mix(PREPARE_GTDB_DB.out.versions)

        // combine contigs and orfs into one channel
        split_metagenome_contigs
            .filter { meta, path ->
                meta.taxon in ['bacteria', 'archaea']
            }
            .map { meta, path -> 
                def cleanMeta = meta.findAll { k,v -> k != 'taxon' }
                [cleanMeta, path]
            }
            .groupTuple(by: 0)
            .combine(
                orfs, by: 0
            )
            .set { contigs_and_orfs_ch }


        EXTRACT_ORFS(contigs_and_orfs_ch)
      
        prot_accession2taxid_ch =Channel.fromPath(file("$baseDir/assets/dummy_file.txt", checkIfExists: true ))
          

        TAXON_SPLIT(
            EXTRACT_ORFS.out.split_orfs,
            PREPARE_GTDB_DB.out.diamond_db,
            PREPARE_GTDB_DB.out.gtdb_taxdump_directory,
            prot_accession2taxid_ch,
            dbtype_ch
        )
        ch_versions = ch_versions.mix(TAXON_SPLIT.out.versions)


        TAXON_SPLIT.out.taxonomically_split_fna.view { meta ->
            println "taxonomically_split_fnabro: ${meta}"
        }

    emit:
        taxonomy                = TAXON_SPLIT.out.taxonomy
        taxonomically_split_fna = TAXON_SPLIT.out.taxonomically_split_fna
        lca                     = TAXON_SPLIT.out.lca
        votes                   = TAXON_SPLIT.out.votes
        taxdump_files           = PREPARE_GTDB_DB.out.gtdb_taxdump_directory
        dbtype                  = dbtype_ch
        versions                = ch_versions

}

