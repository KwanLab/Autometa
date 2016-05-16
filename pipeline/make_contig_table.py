#!/usr/bin/env python

# Program to make a gc, length and coverage table from a spades assembly
# (gets the coverage from the contig name)

import sys
import getopt
from Bio import SeqIO
from Bio.SeqUtils import GC

input_fasta_path = None # -i --input
output_table_path = None # -o --output

opts,args = getopt.getopt(sys.argv[1:],"hi:o",["help", "input=", "output="])

for opt, arg in opts:
	if opt in ('-h', '--help'):
		print 'make_contig_table.py -i <input fasta> -o <output table>'
		sys.exit()
	elif opt in ('-i', '--input'):
		input_fasta_path = arg
	elif opt in ('-o', '--output'):
		output_table_path = arg

output = open( output_table_path, 'w' )

output.write('contig\tlength\tgc\tcov\n')

for seq_record in SeqIO.parse(input_fasta_path, 'fasta'):
	length = len(seq_record)
	gc = GC(str(seq_record.seq))
	contig = str(seq_record.id)
	contigList = contig.split('_')
	cov = contigList[5]
	output.write(contig + '\t' + length + '\t' + gc + '\t' + cov + '\n')

output.close