#!/usr/bin/env python

# Program to make a control assembly fasta for simulated metagenome datasets made by 'make_simulated_metagenome.py'.
# The control fasta contains all the contigs originally obtained from NCBI, and a table is also made that
# tells you which species/assembly each contig came from.  You can align simulated reads to this fasta file in
# order to find out which reads are really unique to certain species.  This information can be used later to determine
# the accuracy of assemblies and binning of the simulated reads.

# For now the paths are hardcoded - this is sort of prototype code

import sys
import subprocess
import pprint

metagenome_table = sys.argv[1]
pp = pprint.PrettyPrinter(indent=4)

# Parse table
base_dir = '/mnt/not_backed_up/ncbi_bacteria_genomes/random_3000_bacteria'
assemblies = {} # Dictionary to hold assembly paths
rows = ((row.rstrip('\n')) for row in open(metagenome_table))
for i,row in enumerate(rows):
	if not i == 0:
		row_list = row.split('\t')
		asm_dir = base_dir + '/' + row_list[0] + '/' + row_list[1]
		# Now we need to find the .fna file
		search_string = asm_dir + '/*.fna'
		ls_output = subprocess.check_output('ls -1 ' + search_string, shell=True)
		fna_list = ls_output.split()
		if not len(fna_list) == 1:
			print 'Error, wrong number of .fna files in ' + asm_dir
			pp.pprint(fna_list)
			sys.exit(2)
		asm_path = fna_list[0]
		#asm_path = asm_dir + '/' + fna_filename
		assemblies[row_list[0]] = asm_path

# Consolidate fastas
# Remember details of contigs
contig_species = {} # Dictionary which will hold the species of each contig
# Output path derived from table filename
table_list = metagenome_table.split('.')
table_list.pop()
output_fasta_filename = ('.').join(table_list) + '_control_contigs.fasta'
output_fasta = open(output_fasta_filename, 'w')
for species in assemblies:
	asm_path = assemblies[species]
	with open(asm_path) as asm:
		for line in asm:
			if line[0] == '>':
				# Get contig name
				name_list = line[1:].split()
				seq_name = name_list[0]
				contig_species[seq_name] = species
				output_fasta.write(line)
			else:
				output_fasta.write(line)
output_fasta.close

# Output contig table
output_table_filename = '.'.join(table_list) + '_control_contigs.tab'
output_table = open(output_table_filename, 'w')
output_table.write('contig\tspecies\n')
for contig in contig_species:
	species = contig_species[contig]
	output_table.write(contig + '\t' + species + '\n')
output_table.close