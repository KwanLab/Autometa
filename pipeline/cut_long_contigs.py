#!/usr/bin/env python

# Program to cut long contigs into uniform size pieces before running vizbin
# USAGE: cut_long_contigs.py <input fasta> <length of pieces in bp> <output fasta>

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pprint
pp = pprint.PrettyPrinter(indent=4)

input_fasta_path = sys.argv[1]
piece_length = sys.argv[2]
output_fasta_path = sys.argv[3]

processed_seq_records = []
for seq_record in SeqIO.parse(input_fasta_path, 'fasta'):
	if len(seq_record) > (3 * piece_length):
		root_id = seq_record.id
		current_seq = seq_record.seq
		section_counter = 0
		for i in range(0, (len(seq_record) - (2 * piece_length)), piece_length):
			section_counter += 1
			if (len(seq_record) - i) > (2 * piece_length):
				new_seq_record = seq_record[i : (i + piece_length - 1)]
				new_seq_record.id = root_id + '_' + str(section_counter)
				new_seq_record.name = root_id + '_' + str(section_counter)
				processed_seq_records.append(new_seq_record)
			else:
				new_seq_record = seq_record[i : (len(seq_record) - 1)]
				new_seq_record.id = root_id + '_' + str(section_counter)
				new_seq_record.name = root_id + '_' + str(section_counter)
				processed_seq_records.append(new_seq_record)

	else:
		processed_seq_records.append(seq_record)

pp.pprint(processed_seq_records)