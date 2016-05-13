#!/usr/bin/env python

# Program to make a simulated lane of illumina reads from a certain number of random bacterial genomes
# The number is expressed as the total length of bacterial genome
# For example - 10,000 Mbp of bacterial genome will be represented in a 400,000,000 paired end 125 bp run
# at a coverage of 10x (all genomes will be of equal coverage)
# The stdout will include read numbers for each component species so that you can subsequently track how 
# these reads are assembled (i.e. to judge the quality of binning procedures)

import sys
import subprocess
from random import randint
import glob
import pprint
import os
import getopt

def fasta_length ( fasta ):
	fasta_length = 0
	lines = (line.rstrip('\n') for line in open(fasta))
	for line in lines:
		if not (line[0] == '>'):
			fasta_length += len(line)
	return fasta_length

pp = pprint.PrettyPrinter(indent=4)

# Get options
length = None
reads = None
bacteria_genome_dir = None
output_prefix = None

#try:
opts,args = getopt.getopt(sys.argv[1:],"hl:r:i:o:",["help", "length=", "reads=", "genomedir=", "outprefix="])
#except getopt.GetoptError:
#	print 'make_simulated_metagenome.py -l <total bacterial genome length in Mbp to include> -r <total reads to output> -i <bacterial genomes dir> -o <output prefix, including path>'
#	sys.exit(2)
for opt, arg in opts:
	if opt in ('-h', '--help'):
		print 'make_simulated_metagenome.py -l <total bacterial genome length in Mbp to include> -r <total reads to output> -i <bacterial genomes dir> -o <output prefix, including path>'
		sys.exit()
	elif opt in ('-l', '--length'):
		length = arg
	elif opt in ('-r', '--reads'):
		reads = arg
	elif opt in ('-i', '--genomedir'):
		bacteria_genome_dir = arg
	elif opt in ('-o', '--outputprefix'):
		output_prefix = arg

print 'length: ' + length
print 'reads: ' + reads
print 'genome dir: ' + bacteria_genome_dir
print 'output prefix: ' + output_prefix

# Fix types
length = float(length)
reads = int(reads)

# Work out coverage level for art_illumina
# C = [(LN)/G]/2 
# C = coverage
# L = read length (total of paired reads)
# G = genome size in bp
coverage = ((250 * reads) / (length * 1000000))

# First get list of input bacterial species
#bacteria_genome_dir = '/mnt/not_backed_up/ncbi_bacteria_genomes/random_3000_bacteria'
ls_output = subprocess.check_output(['ls', '-1', bacteria_genome_dir])
bacteria_list = ls_output.split()
highest_index = len(bacteria_list) - 1

# Now make a random list of bacteria, whose genomes add up to ~length
total_length = 0
random_bacteria = dict() # Will store dictionaries of assembly path and read range
random_bacteria_indices = dict() # Will store indexes of the bacteria picked
print 'Picking bacteria...'
while (total_length < (length * 1000000)):
	i = randint(0, highest_index)
	while (i in random_bacteria_indices):
		i = randint(0, highest_index)

	random_bacteria_indices[i] = 1

	# First get list of dirs within the bacteria dir - should be only one there
	bacterium_dir = bacteria_genome_dir + '/' + bacteria_list[i]
	bacterium_ls_output = subprocess.check_output(['ls', '-1', bacterium_dir])
	bacterium_ls_list = bacterium_ls_output.split()
	if not (len(bacterium_ls_list) == 1):
		print 'There is more than one dir/file in ' + bacterium_dir
		sys.exit()
	asm_dir = bacterium_dir + '/' + bacterium_ls_list[0]
	print asm_dir

	# Now we need to work out if the assembly has already been uncompressed, and uncompress it if necessary
	os.chdir(asm_dir)
	compressed_list = glob.glob('*_genomic.fna.gz')
	uncompressed_asm = ''
	if (len(compressed_list) == 1):
		# Uncompress the file
		subprocess.call(['gunzip', compressed_list[0]])
		# Get uncompressed filename
		filename_list = compressed_list[0].split('.')
		filename_list.pop()
		uncompressed_asm = '.'.join(filename_list)
	elif (len(compressed_list) == 0):
		# In this case there will hopefully already be an uncompressed file
		uncompressed_list = glob.glob('*_genomic.fna')
		if (len(uncompressed_list) == 1):
			uncompressed_asm = uncompressed_list[0]
		else:
			print asm_dir + ": can't find *_genomic.fna.gz or *_genomic.fna"
			sys.exit()
	else:
		print asm_dir + ': there seems to be more than one *_genomic.fna.gz files in this directory'
		sys.exit()

	# Now we need to calculate the assembly's total length
	asm_path = asm_dir + '/' + uncompressed_asm
	asm_length = fasta_length(asm_path)
	total_length += asm_length

	# Add the assembly to the dictionary
	tempDict = dict()
	tempDict['asm_path'] = asm_path
	random_bacteria[bacteria_list[i]] = tempDict
	random_bacteria[bacteria_list[i]]['length'] = asm_length

#pp.pprint(random_bacteria)
# Print out a table for debugging purposes
print 'bacterium\tasm\tlength'
for bacterium in random_bacteria:
	asm = random_bacteria[bacterium]['asm_path']
	length = random_bacteria[bacterium]['length']
	print bacterium + '\t' + asm + '\t' + str(length)
#sys.exit()


# Now we need to run art_illumina for all genomes and create a simulated lane of illumina reads
R1_output_path = output_prefix + '_R1.fastq'
R2_output_path = output_prefix + '_R2.fastq'
R1 = open(R1_output_path, 'w')
R2 = open(R2_output_path, 'w')
readcounter = 1

for bacterium in random_bacteria:
	asm_path = random_bacteria[bacterium]['asm_path']
	asm_list = asm_path.split('/')
	asm_list.pop()
	asm_dir = '/'.join(asm_list)

	os.chdir(asm_dir)
	print 'Simulating reads for' + asm_dir
	subprocess.call(['art_illumina', '-p', '-ss', 'HS25', '-l', '125', '-f', str(coverage), '-o', 'simulated_reads', '-m', '275', '-s', '90', '-i', asm_path])

	reads1_path = asm_dir + '/simulated_reads1.fq'
	reads2_path = asm_dir + '/simulated_reads2.fq'

	R2_read_number = readcounter
	start_read = readcounter

	# First do R1
	print 'Writing to R1'
	line_counter = 1
	with open(reads1_path) as R1_input:
		for line in R1_input:
			if line_counter == 1:
				R1.write('@read_' + str(readcounter) + '/1\n')
				readcounter += 1
			else:
				R1.write(line)
			line_counter += 1

			if line_counter == 5:
				line_counter = 1

	end_read = readcounter - 1

	# Now do R2
	print 'Writing to R2'
	line_counter = 1
	with open(reads2_path) as R2_input:
		for line in R2_input:
			if line_counter == 1:
				R2.write('@read_' + str(R2_read_number) + '/2\n')
				R2_read_number += 1
			else:
				R2.write(line)
			line_counter += 1

			if line_counter == 5:
				line_counter = 1

	# Record read range
	random_bacteria[bacterium]['start_read'] = start_read
	random_bacteria[bacterium]['end_read'] = end_read

	# We need to now clean up the created files, to save disk space
	subprocess.call(['rm', 'simulated_reads1.fq'])
	subprocess.call(['rm', 'simulated_reads2.fq'])
	subprocess.call(['rm', 'simulated_reads1.aln'])
	subprocess.call(['rm', 'simulated_reads2.aln'])

R1.close
R2.close

# Now output table of read ranges
table_output_path = output_prefix + '.tab'
table = open(table_output_path, 'w')
table.write('bacterium\tassembly\tstart_read\tend_read\n')
for bacterium in random_bacteria:
	asm_path = random_bacteria[bacterium]['asm_path']
	asm_path_list = asm_path.split('/')
	asm_name = asm_path_list[-2]
	start_read = random_bacteria[bacterium]['start_read']
	end_read = random_bacteria[bacterium]['end_read']

	table.write(bacterium + '\t' + asm_name + '\t' + str(start_read) + '\t' + str(end_read) + '\n')

table.close