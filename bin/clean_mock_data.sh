#!/bin/bash

zcat **/*.fna.gz > combined_nucleotide.fna

splitter -size 5000 combined_nucleotide.fna -outseq metagenome.fna

gzip metagenome.fna
gzip combined_nucleotide.fna

# get assembly accessions/locus
zgrep ">" **/*.fna.gz  | sed 's/\s.*$//' | sed -E 's/[^\/]*.>//' > assembly_to_locus.txt #TODO: delimit better
# get assembly accessions
zgrep ">" **/*.fna.gz  | cut -d_ -f1,2  | uniq > assemblies.txt

# Add fake spades coverage to deflines
zcat "metagenome.fna.gz" | awk '/^>/ {$0=$1} 1' | sed 's/>.*/&_length_1_cov_1/' | gzip > fake_spades.fna.gz


