#!/usr/bin/env python

# Program to make a gc, length and coverage table from a spades assembly
# (gets the coverage from the contig name)

import sys
import getopt
from Bio import SeqIO
from Bio.SeqUtils import GC

input_fasta_path = sys.argv[1]
output_table_path = sys.argv[2]

output = open( output_table_path, 'w' )

output.write('contig\tlength\tgc\tcov\n')

for seq_record in SeqIO.parse(input_fasta_path, 'fasta'):
	length = str(len(seq_record))
	gc = str(GC(str(seq_record.seq)))
	contig = str(seq_record.id)
	contigList = contig.split('_')
	cov = str(contigList[5])
	output.write(contig + '\t' + length + '\t' + gc + '\t' + cov + '\n')

output.close