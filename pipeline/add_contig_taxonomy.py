#!/usr/bin/env python

# Program that adds contig taxonomy information to a table made by make_contig_table.py
# Uses protein taxonomy information to construct contig taxonomy
# Algorithm:
# In descending order of species, genus, family, order, class, phylum, superkingdom:
#    In descending order of votes:
#        If classification shares common ancestry with majority of other proteins, accept result
#    If no result, move up to next taxonomic level

import sys
from time import gmtime, strftime
import pprint
pp = pprint.PrettyPrinter(indent=4)

tax_table_path = sys.argv[1]
taxdump_dir_path = sys.argv[2]
output_file_path = sys.argv[3]

# Process NCBI taxdump files
names_dmp_path = taxdump_dir_path + '/names.dmp'
nodes_dmp_path = taxdump_dir_path + '/nodes.dmp'

taxids = {}
print strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid names'
with open(names_dmp_path) as names_dmp:
	for line in names_dmp:
		line_list = line.rstrip('\n').split('|')
		# Remove trailing and leading spaces
		for i,value in line_list:
			line_list[i] = value.strip

		# Only add scientific name entries
		if line_list[3] == 'scientific name':
			line_list[1].replace(' ', '_')
			taxids[line_list[0]] = { 'name': line_list[1] }

print strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid nodes'
with open(nodes_dmp_path) as nodes_dmp:
	for line in nodes_dmp:
		line_list = line.rstrip('\n').split('|')
		# Remove trailing and leading spaces
		for i,value in line_list:
			line_list[i] = value.rstrip

		taxids[ line_list[0] ][ 'parent' ] = line_list[1]
		taxids[ line_list[0] ][ 'rank' ] = line_list[2]

pp.pprint(taxids)
