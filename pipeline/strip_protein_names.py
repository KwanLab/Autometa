#!/usr/bin/env python

# Program to de-identify protein names in a blast table, so that all proteins from the same contig are only identified by the contig name
# Why do we want to do this?  So that the data can be passed through a lowest common ancestor calculator (e.g. blast2lca) to obtain
# taxonomy of all contigs
# Note - this program is very dumb and doesn't do a whole lot of checking

import sys

input_table_path = sys.argv[1]
output_table_path = sys.argv[2]

output_table = open(output_table_path, 'w')

with open(input_table_path) as input_table:
	for line in input_table:
		line_list = line.split('\t')
		query_name = line_list[0]
		line_list.pop(0)
		query_list = query_name.split('_')
		query_list.pop()
		new_query_name = ('_').join(query_list)
		new_line = new_query_name + '\t' + ('\t').join(line_list)
		output_table.write(new_line)

output_table.close