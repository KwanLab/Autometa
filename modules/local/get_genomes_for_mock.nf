process GET_GENOMES_FOR_MOCK {
    def genome_count = task.ext.args2.tokenize('|').size()
    tag "fetching ${genome_count} genomes"

    cache 'lenient'

    conda "bioconda::emboss=6.6.0"
    container "jasonkwan/autometa-nf-modules-get_genomes_for_mock:${params.autometa_image_tag}"

    output:
        path "metagenome.fna.gz", emit: metagenome
        path "combined_nucleotide.fna.gz", emit: combined_nucleotide
        path "fake_spades.fna.gz", emit: fake_spades_coverage
        path "assembly_to_locus.txt", emit: assembly_to_locus
        path "assemblies.txt", emit: assemblies
        path "assembly_report.txt", emit: assembly_report

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
    """
    curl -s ${args} > assembly_report.txt

    grep -E ${args2} assembly_report.txt |\\
        awk -F '\\t' '{print \$20}' |\\
        sed 's,https://,rsync://,' |\\
            xargs -n 1 -I {} \
                rsync -am \
                    --exclude='*_rna_from_genomic.fna.gz' \
                    --exclude='*_cds_from_genomic.fna.gz' \
                    --include="*_genomic.fna.gz" \
                    --include="*_protein.faa.gz" \
                    --include='*/' \
                    --exclude='*' {} .

    # "clean_mock_data.sh" is here: ~/Autometa/bin/clean_mock_data.sh
    clean_mock_data.sh
    """
}
