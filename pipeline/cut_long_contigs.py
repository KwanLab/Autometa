#!/usr/bin/env python

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

# Program to cut long contigs into uniform size pieces before running vizbin
# USAGE: cut_long_contigs.py <input fasta> <length of pieces in bp> <output fasta>

# Note: there is currently a bug in this program that causes a little bit of sequence to be lost every time it is run.
# I'm not sure why this happens but we should try to work out what is going on.

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pprint
pp = pprint.PrettyPrinter(indent=4)

input_fasta_path = sys.argv[1]
piece_length = int(sys.argv[2])
output_fasta_path = sys.argv[3]

processed_seq_records = []
for seq_record in SeqIO.parse(input_fasta_path, 'fasta'):
	print('ID: ' + seq_record.id)
	print('Length: ' + str(len(seq_record)))
	if len(seq_record) > (3 * piece_length):
		root_id = seq_record.id
		current_seq = seq_record.seq
		section_counter = 0
		for i in range(1, len(seq_record), piece_length):
			section_counter += 1
			if (len(seq_record) - i) > (2 * piece_length):
				print ('Cut start: ' + str(i) + ' Cut end: ' + str(i + piece_length - 1))
				new_seq_record = seq_record[i : (i + piece_length - 1)]
				new_seq_record.id = root_id + '_' + str(section_counter)
				new_seq_record.name = root_id + '_' + str(section_counter)
				processed_seq_records.append(new_seq_record)
			else:
				print ('Cut start: ' + str(i) + ' Cut end: ' + str(len(seq_record)))
				new_seq_record = seq_record[i : len(seq_record)]
				new_seq_record.id = root_id + '_' + str(section_counter)
				new_seq_record.name = root_id + '_' + str(section_counter)
				processed_seq_records.append(new_seq_record)
				break
	else:
		processed_seq_records.append(seq_record)

#pp.pprint(processed_seq_records)
print 
print ('Processed records: ')
print
for seq_record in processed_seq_records:
	print ('ID: ' + seq_record.id)
	print ('Length: ' + str(len(seq_record)))


# Now write the seq_records to a fasta file
SeqIO.write(processed_seq_records, output_fasta_path, 'fasta')