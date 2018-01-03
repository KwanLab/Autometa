#!/usr/bin/env python

# Program to make a gc, length and coverage table from an assembly
# If you do not specify a coverage table, it attempts to use the contig name

import argparse
import os
import sys
import getopt
from Bio import SeqIO
from Bio.SeqUtils import GC

#argument parser
parser = argparse.ArgumentParser(description="Script to make a gc, length and coverage table from an assembly",\
  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a', '--assembly', metavar='<assembly.fasta>', help='Path to assembly fasta', required=True)
parser.add_argument('-c', '--coverage', metavar='<coverage.tab>', help='Path to coverage table made by calculate_read_coverage.py \
	 [if not supplied coverage will be inferred from contig names]')
parser.add_argument('-o', '--output', metavar='<output.tab>', help='Path to create output table', required=True)

args = vars(parser.parse_args())

fasta_assembly_path = os.path.abspath(args['a'])
coverage_table_path = args['c']
output_table_path = args['o']

use_coverage_table = False

if coverage_table_path is not None:
	use_coverage_table = True
	# We need to parse the coverage table now
	coverages = dict() # Will be keyed by contig name
	if os.path.isfile(coverage_table_path):
		with open(coverage_table_path) as table:
			for i,line in enumerate(table):
				if i > 0:
					line_list = line.rstrip().split('\t')
					contig = line_list[0]
					coverage = float(line_list[1])
					coverages[contig] = coverage
	else:
		print('Error, couldnt find file ' + coverage_table_path + ' or it is unreadable')
		quit()

#input_fasta_path = sys.argv[1]
#output_table_path = sys.argv[2]

output = open( output_table_path, 'w' )

output.write('contig\tlength\tgc\tcov\n')

for seq_record in SeqIO.parse(input_fasta_path, 'fasta'):
	length = str(len(seq_record))
	gc = str(GC(str(seq_record.seq)))
	contig = str(seq_record.id)
	if use_coverage_table:
		if contig in coverages:
			cov = coverages[contig]
		else:
			print('Error, ' + contig + ' not found in the coverage table')
			quit()
	else:
		# Do a format check to make sure the contig name is right
		contigList = contig.split('_')
		if len(contigList) == 6 and contigList[0] == 'NODE' and contigList[2] == 'length' and contigList[4] == 'cov':
			cov = str(contigList[5])
		else:
			print('Error, ' + contig + ' not the right format to extract coverage from sequence name')
			quit()
	output.write(contig + '\t' + length + '\t' + gc + '\t' + cov + '\n')

output.close