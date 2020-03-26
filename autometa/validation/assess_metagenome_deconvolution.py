#!/usr/bin/env python2.7

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

# Program to assess the quality of deconvolution of a simulated or defined metagenome dataset
# Inputs:
# 	1. Sam file of reads aligned to the known reference genomes
#	2. Table of species classifications of each contig in the reference genome fasta
#	3. Table of read ranges for each species
#	4. Sam file of reads aligned to the deconvoluted assembly
#	5. Table of bin classifications of each contig in the assembly

# Outputs:
#	1. 'Binning accuracy' table, header: bin\tgenome\tpercent
#		Tells you the percent of each bin that belongs to different genomes
#	2. 'Binning recovery' table, header: genome\tbin\tpercent
#		Tells you the percent of each genome found in a bin
#	3. 'Chimera' table, header: contig\tgenome\tpercent
#		Tells you the percent of each contig that belongs to different genomes

import sys
import argparse 
import gzip
from time import *
from tqdm import *
import subprocess
import pprint
import os.path
pp = pprint.PrettyPrinter(indent=4)

def is_alignment_congruent_with_ref(read_name, contig_aligned_to, read_ranges, contig_species):
	species = contig_species[contig_aligned_to]
	if does_read_belong(read_name, species, read_ranges):
		return True
	else:
		return False

def does_read_belong(read_name, species, read_ranges):
	# Read name assumed to be in the form read_235
	read_list = read_name.split('_')
	read_number = read_list[1]
	species_start_read = read_ranges[species]['start']
	species_end_read = read_ranges[species]['end']
	if read_number >= species_start_read and read_number <= species_end_read:
			return True
	else:
		return False

def species_classification_for_read(read_name, read_ranges, excluded_reads):
	classified_species = None
	for species in read_ranges:
		if does_read_belong(read_name, species, read_ranges):
			if read_name not in excluded_reads:
				classified_species = species
	return classified_species

def get_species_percents(read_counts):
	percents = {}
	total_read_count = 0
	for species in read_counts:
		total_read_count += read_counts[species]
	if total_read_count == 0:
		# This can happen if there is missing data (e.g. a junk bin)
		return None
	for species in read_counts:
		percent = (float(read_counts[species])/float(total_read_count))*100
		percents[species] = percent
	return percents

parser = argparse.ArgumentParser(description='Script to assess the metagenome deconvolution method')
parser.add_argument('-r','--refsam', help='sam alignment to ref genomes', required=True)
parser.add_argument('-s', '--refspecies', help='ref contig species table', required=True)
parser.add_argument('-q', '--readranges', help='ref read table', required=True)
parser.add_argument('-a', '--asmsam', help='sam alignment to assembly', required=True)
parser.add_argument('-b', '--bintable', help='bin classification table for assembly', nargs='+', required=True)
parser.add_argument('-c', '--column', help='bin column name', default = 'db.cluster')
parser.add_argument('-o', '--outputprefix', help='bin column name', required=True)
args = vars(parser.parse_args())

ref_sam_path = args['refsam']
ref_species_table_path = args['refspecies']
ref_read_ranges_table_path = args['readranges']
asm_sam_path = args['asmsam']
bin_classifications_table_paths =  args['bintable']
output_prefix = args['outputprefix']
bin_column = args['column']

if not os.path.exists(ref_sam_path):
	print 'Error, could not find ' + ref_sam_path
	sys.exit(2)
print 'Reference SAM: ' + ref_sam_path
if not os.path.exists(ref_species_table_path):
	print 'Error, could not find ' + ref_species_table_path
	sys.exit(2)
print 'Contig species table: ' + ref_species_table_path
if not os.path.exists(ref_read_ranges_table_path):
	print 'Error, could not find ' + ref_read_ranges_table_path
	sys.exit(2)
print 'Read ranges table: ' + ref_read_ranges_table_path
if not os.path.exists(asm_sam_path):
	print 'Error, could not find ' + asm_sam_path
	sys.exit(2)
print 'Assembly SAM: ' + asm_sam_path
for bin_classifications_table_path in bin_classifications_table_paths:
	if not os.path.exists(bin_classifications_table_path):
		print 'Error, could not find ' + bin_classifications_table_path
		sys.exit(2)
print 'Bin classifications table(s): ' + (' ').join(bin_classifications_table_paths)
# To do - move column detection here so as not to waste time if the tables are not formatted right
print 'Bin column: ' + bin_column
output_list = output_prefix.split('/')
output_list.pop()
output_dir = ('/').join(output_list)
if not os.path.isdir(output_dir):
	print 'Error, ' + output_dir + ' is not a valid directory for output'
	sys.exit(2)
print 'Output prefix: ' + output_prefix
print strftime("%Y-%m-%d %H:%M:%S")
print 

# 1. Parse read ranges table, so that we can spot non-unique reads in the reference alignment
print strftime("%Y-%m-%d %H:%M:%S") + ' Parsing read ranges table...'
wc_output = subprocess.check_output(['wc', '-l', ref_read_ranges_table_path])
wc_list = wc_output.split()
number_of_lines = int(wc_list[0])

ranges = {} # Dictionary of dictionaries, keyed by species
range_table_rows = ((row.rstrip('\n')) for row in open(ref_read_ranges_table_path))
last_read = 0
for i,row in enumerate(tqdm(range_table_rows, total=number_of_lines)):
	if not i == 0:
		row_list = row.split('\t')
		ranges[row_list[0]] = { 'start':row_list[2], 'end':row_list[3] }
		last_read = row_list[3]

# 2. Go through reference contig table, and remember which species each contig belongs to
print strftime("%Y-%m-%d %H:%M:%S") + ' Parsing contig species table...'
wc_output = subprocess.check_output(['wc', '-l', ref_species_table_path])
wc_list = wc_output.split()
number_of_lines = int(wc_list[0])
number_of_contigs = number_of_lines

number_sam_lines = (int(last_read)*2) + number_of_lines # We assume paired-end reads here

species = {} # Dictionary, keyed by contig, stores species
species_table_rows = ((row.rstrip('\n')) for row in open(ref_species_table_path))
for i,row in enumerate(tqdm(species_table_rows, total=number_of_lines)):
	if not i == 0:
		row_list = row.split('\t')
		species[row_list[0]] = row_list[1]

# 3. Go through reference sam file, and flag non-unique reads
# Non-unique reads are spotted in the sam file as reads that occur in genomes other than their originating genome
# (because they must be in their originating genome, therefore seeing one outside means they occur in at least two)
# Note: this does not flag up reads that occur more than once in their originating genomes
print strftime("%Y-%m-%d %H:%M:%S") + ' Finding non-unique reads in reference SAM...'
non_unique_reads = {} # Dictionary that just contains reads found in more than one genome

# If sam file is a gz file, use gzip, otherwise normal open
if ref_sam_path[-3:] == '.gz':
	with gzip.open(ref_sam_path, 'rb') as ref_sam:
		for line in tqdm(ref_sam, total=number_sam_lines):
			# Skip header section, where lines begin with '@'
			if not line[0] == '@':
				line_list = line.split('\t')
				read_name = line_list[0]
				contig_name = line_list[2]
				contig_species = species[contig_name]
				# Find if read "belongs" to the species that this contig is from
				if not does_read_belong(read_name, contig_species, ranges):
					non_unique_reads[read_name] = 1
else:
	with open(ref_sam_path) as ref_sam:
		for line in tqdm(ref_sam, total=number_sam_lines):
			# Skip header section, where lines begin with '@'
			if not line[0] == '@':
				line_list = line.split('\t')
				read_name = line_list[0]
				contig_name = line_list[2]
				contig_species = species[contig_name]
				# Find if read "belongs" to the species that this contig is from
				if not does_read_belong(read_name, contig_species, ranges):
					non_unique_reads[read_name] = 1

# 3a. Make a data structure that records the number of unique reads for each bin.
print strftime("%Y-%m-%d %H:%M:%S") + ' Working out how many unique reads there are per genome...'
number_of_unique_reads = {} # Dictionary keyed by bin name
number_of_genomes = len(ranges)

for genome in tqdm(ranges, total=number_of_genomes):
	start_read = int(ranges[genome]['start'])
	end_read = int(ranges[genome]['end'])
	counter = start_read
	unique_reads = 0
	while counter <= end_read:
		read_name = 'read_' + str(counter)
		if read_name not in non_unique_reads:
			unique_reads += 1
		counter += 1
	number_of_unique_reads[genome] = unique_reads

# 4. Go through assembly sam file, and count read classifications for each contig
print strftime("%Y-%m-%d %H:%M:%S") + ' Parsing assembly SAM, counting species reads...'
contig_classifications = {} # Dictionary of dictionaries, which will hold running tallies of reads assigned to different species

# If sam file is a gz file, use gzip
if asm_sam_path[-3:] == '.gz':
	with gzip.open(asm_sam_path, 'rb') as asm_sam:
		for line in tqdm(asm_sam, total=number_sam_lines):
			# Skip header section, where lines begin with '@'
			if not line[0] == '@':
				line_list = line.split('\t')
				read_name = line_list[0]
				contig_name = line_list[2]
				read_species = species_classification_for_read(read_name, ranges, non_unique_reads) # Returns None if a non-unique read
				if read_species:
					if contig_name in contig_classifications:
						if read_species in contig_classifications[contig_name]:
							contig_classifications[contig_name][read_species] += 1
						else:
							contig_classifications[contig_name][read_species] = 1
					else:
						contig_classifications[contig_name] = { read_species: 1 }
else:
	with open(asm_sam_path) as asm_sam:
		for line in tqdm(asm_sam, total=number_sam_lines):
			# Skip header section, where lines begin with '@'
			if not line[0] == '@':
				line_list = line.split('\t')
				read_name = line_list[0]
				contig_name = line_list[2]
				read_species = species_classification_for_read(read_name, ranges, non_unique_reads) # Returns None if a non-unique read
				if read_species:
					if contig_name in contig_classifications:
						if read_species in contig_classifications[contig_name]:
							contig_classifications[contig_name][read_species] += 1
						else:
							contig_classifications[contig_name][read_species] = 1
					else:
						contig_classifications[contig_name] = { read_species: 1 }

# 5. We now have enough information to write a table showing how chimeric contigs are
# Output table in the format contig\tgenome\treads\tpercent
chimera_table_path = output_prefix + '_chimera_table'
print strftime("%Y-%m-%d %H:%M:%S") + ' Writing chimera table ' + chimera_table_path + '...'
chimera_table = open(chimera_table_path, 'w')
chimera_table.write('contig\tgenome\treads\tpercent\n')
for contig in tqdm(contig_classifications, total=number_of_contigs):
	percents = get_species_percents(contig_classifications[contig])
	for species in contig_classifications[contig]:
		chimera_table.write(contig + '\t' + species + '\t' + str(contig_classifications[contig][species]) + '\t' + str(percents[species]) + '\n')
chimera_table.close

# 5a. Make a table that determines the yield of genome reads that align to contigs (i.e. judge the quality of the assembly)
aligned_unique_reads = {} # Dictionary that will store totals of aligned reads for each species
reads_in_pure_contigs = {} # Dictionary that will store totals of aligned reads in 'pure' > 90% contigs
reads_in_chimeric_contigs = {} # Dictionary that will store totals of aligned reads in 'impure' < 90% contigs
# Note, according to the above scheme, the same contig can be considered pure and impure from the point of view of different genomes
# i.e. - a contigs that is 95% E. coli and 5% S. aureus would be a pure E. coli contig and an impure S. aureus contig

asm_recovery_table_path = output_prefix + '_asm_recovery_table'
print strftime("%Y-%m-%d %H:%M:%S") + ' Writing assembly recovery table ' + asm_recovery_table_path + '...'
for contig in contig_classifications:
	percents = get_species_percents(contig_classifications[contig])
	for species in percents:
		if species not in aligned_unique_reads:
			aligned_unique_reads[species] = 0
		if species not in reads_in_pure_contigs:
			reads_in_pure_contigs[species] = 0
		if species not in reads_in_chimeric_contigs:
			reads_in_chimeric_contigs[species] = 0

		aligned_unique_reads[species] += contig_classifications[contig][species]

		if percents[species] >= 90:
			reads_in_pure_contigs[species] += contig_classifications[contig][species]
		else:
			reads_in_chimeric_contigs[species] += contig_classifications[contig][species]

# Write out table
asm_recovery_table = open(asm_recovery_table_path, 'w')
asm_recovery_table.write('species\taligned_unique_reads\treads_in_pure_contigs\treads_in_chimeric_contigs\n')
for species in aligned_unique_reads:
	asm_recovery_table.write(species + '\t' + str(aligned_unique_reads[species]) + '\t' + str(reads_in_pure_contigs[species]) + '\t' + str(reads_in_chimeric_contigs[species]) + '\n')
asm_recovery_table.close

# As we iterate through the supplied tables, we will claculate the 'clustering quotient' for each one, and hold it in a data structure
clustering_quotients = {} # Keyed by bin_classification_table_path
bin_accuracy_averages = {}
genome_recovery_averages = {}
# Total_bin_genome_recovery - this is a measure of the effectiveness of the deconvolution process
# Algorithm
# For each table
#     For each genome
#         Determine which bin contains the highest percentage of this genome, as long as it is not claimed
#         Add this percentage to a total, mark the bin as claimed
# The total is stored in total_bin_genome_recovery, under the table path
total_bin_genome_recovery = {} 

##### Here we iterate over all entries in bin_classifications_table_paths
for bin_classifications_table_path in bin_classifications_table_paths:
	print 'Considering ' + bin_classifications_table_path
	pathList = bin_classifications_table_path.split('/')
	#filename = None
	#if '.' in pathList[-1]:
	#	filenameList = pathList[-1].split('.')
	#	filenameList.pop()
	#	filename = ('.').join(filenameList)
	#else: 
	filename = pathList[-1]

	current_output_prefix = output_prefix + '_' + filename

	# 6. We need to go through the bin table to make a datastructure containing the classification of each contig
	print strftime("%Y-%m-%d %H:%M:%S") + ' Making bin datastructure...'
	contig_bins = {} # Dictionary, keyed by contig, stores bin classifications
	bin_table_rows = None
	with open(bin_classifications_table_path) as bin_classifications_table:
		bin_table_rows = bin_classifications_table.read().splitlines()

	# Find out which column to use for bin classification 
	bin_column_index = None
	number_found = 0
	first_line_list = bin_table_rows[0].split('\t')
	for i,value in enumerate(first_line_list):
		if value == bin_column:
			bin_column_index = i
			number_found += 1
	if number_found > 1:
		print 'Error, bin table has more than one column headed ' + bin_column
		sys.exit(2)
	if not bin_column_index:
		print 'Error, could not find column ' + bin_column + ' in bin table ' + bin_classifications_table_path
		sys.exit(2)

	contig_column_index = None
	number_found = 0
	for i,value in enumerate(first_line_list):
		if value == 'contig':
			contig_column_index = i
			number_found += 1
	if number_found > 1:
		print 'Error, bin table has more than one contig column'
		sys.exit(2)
	if contig_column_index is None:
		print 'Error, could not find contig column in bin table'
		sys.exit(2)

	for i,row in enumerate(bin_table_rows):
		if i != 0:
			row_list = row.split('\t')
			contig = row_list[contig_column_index]
			bin_name = row_list[bin_column_index]
			contig_bins[contig] = bin_name

	# Now make a data structure that totals up the species reads for each bin
	bin_classifications = {} # Dictionary keyed by bins, then species
	total_contigs = len(contig_bins.keys())
	for contig in tqdm(contig_bins, total=total_contigs):
		current_bin = contig_bins.get(contig)
		if current_bin is None:
			print 'Contig ' + contig + ' not found in contig_bins'
			pp.pprint(contig_bins)
		if current_bin not in bin_classifications:
			bin_classifications[current_bin] = {}

		# Sometimes contigs that are in the bin_classifications table do not exist in the alignments
		# - perhaps this means these are junk contigs?
		# Anyway, we have to check here that these exist

		if contig in contig_classifications:
			for species in contig_classifications[contig]:
				number_reads = contig_classifications[contig][species]
				if species in bin_classifications[current_bin]:
					bin_classifications[current_bin][species] += number_reads
				else:
					bin_classifications[current_bin][species] = number_reads

	# 7. Make 'Binning accuracy' table, header: bin\tgenome\treads\tpercent
	bin_accuracy_values = []
	bin_accuracy_table_path = current_output_prefix + '_bin_accuracy_table'
	print strftime("%Y-%m-%d %H:%M:%S") + ' Writing binning accuracy table ' + bin_accuracy_table_path + '...'
	bin_accuracy_table = open(bin_accuracy_table_path, 'w')
	bin_accuracy_table.write('bin\tgenome\treads\tpercent\n')
	for bin_name in bin_classifications:
		percents = get_species_percents(bin_classifications[bin_name])

		if percents is not None:
			percent_list = percents.values()
			percent_list.sort(reverse=True)
			bin_accuracy_values.append(percent_list[0])
			for species in bin_classifications[bin_name]:
				bin_accuracy_table.write(bin_name + '\t' + species + '\t' + str(bin_classifications[bin_name][species]) + '\t' + str(percents[species]) + '\n')
	bin_accuracy_table.close

	# 8. Make a 'Binning recovery' table, header: genome\tbin\treads\tpercent
	print 'Counting genome reads in bins...'
	genome_recovery_values = []
	genome_reads_in_bins = {}
	for bin_name in bin_classifications:
		for species in bin_classifications[bin_name]:
			if species not in genome_reads_in_bins:
				genome_reads_in_bins[species] = {}

			if bin_name in genome_reads_in_bins[species]:
				genome_reads_in_bins[species][bin_name] += bin_classifications[bin_name][species]
			else:
				genome_reads_in_bins[species][bin_name] = bin_classifications[bin_name][species]

	bin_recovery_table_path = current_output_prefix + '_bin_recovery_table'
	print strftime("%Y-%m-%d %H:%M:%S") + ' Writing binning recovery table ' + bin_recovery_table_path + '...'
	bin_recovery_table = open(bin_recovery_table_path, 'w')
	bin_recovery_table.write('genome\tbin\treads\tpercent_of_input\tpercent_of_aligned_to_asm\n')
	genome_recovery_total = 0
	claimed_bins = {}
	for species in genome_reads_in_bins:
		# Percents are calculated based on the previously calculated number of unique reads per reference genome
		# because - not all the genome might end up in bins/assembled
		percents_unique_input = {} # Percent expressed as a fraction of the unique reads aligned to reference genomes
		percents_aligned_to_contigs = {} # Percent expressed as a fraction of the unique reads aligned to the assembly
		for bin_name in genome_reads_in_bins[species]:
			number_of_reads = genome_reads_in_bins[species][bin_name]
			percent_unique_input = (float(number_of_reads) / float(number_of_unique_reads[species]))*50 # Here we assume paired reads
			percent_aligned_to_contigs = (float(number_of_reads) / float(aligned_unique_reads[species]))*50
			percents_unique_input[bin_name] = percent_unique_input
			percents_aligned_to_contigs[bin_name] = percent_aligned_to_contigs

		recovery_values_list = percents_aligned_to_contigs.values()
		recovery_values_list.sort(reverse=True)
		genome_recovery_values.append(recovery_values_list[0])

		# Add to the genome recovery total
		sorted_bins = sorted(percents_aligned_to_contigs, key=percents_aligned_to_contigs.get, reverse=True)
		for bin_name in sorted_bins:
			if bin_name not in claimed_bins:
				genome_recovery_total += percents_aligned_to_contigs[bin_name]
				claimed_bins[bin_name] = 1
				break

		for bin_name in genome_reads_in_bins[species]:
			bin_recovery_table.write(species + '\t' + bin_name + '\t' + str(genome_reads_in_bins[species][bin_name]) + '\t' + str(percents_unique_input[bin_name]) + '\t' + str(percents_aligned_to_contigs[bin_name]) + '\n')
	bin_recovery_table.close

	total_bin_genome_recovery[bin_classifications_table_path] = genome_recovery_total

	# Work out clustering quotient

	# First we have to pad out genome_recovery_values if some of the original species are not missing
	while len(genome_recovery_values) < len(ranges.keys()):
		genome_recovery_values.append(float(0))

	bin_accuracy_average = float(sum(bin_accuracy_values)) / len(bin_accuracy_values)
	genome_recovery_average = float(sum(genome_recovery_values)) / len(genome_recovery_values)

	cluster_quotient = bin_accuracy_average - genome_recovery_average
	clustering_quotients[bin_classifications_table_path] = cluster_quotient
	bin_accuracy_averages[bin_classifications_table_path] = bin_accuracy_average
	genome_recovery_averages[bin_classifications_table_path] = genome_recovery_average


# Print out table of cluster quotients
clustering_table_path = output_prefix + '_clustering_quotients'
clustering_table = open(clustering_table_path, 'w')
clustering_table.write('table\tav_bin_accuracy\tav_genome_recovery\trecovery_total\n')
for table in clustering_quotients:
	clustering_table.write(table + '\t' + str(bin_accuracy_averages[table]) + '\t' + str(genome_recovery_averages[table]) + '\t' + str(total_bin_genome_recovery[table]) + '\n')

clustering_table.close
