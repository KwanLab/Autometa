#!/usr/bin/env python

import sys
import argparse 

parser = argparse.ArgumentParser(description='Script to pick the best dbscan table using a table made by assess_clustering.py')
parser.add_argument('-i','--inputtable', help='HMM table created by make_marker_table.py', required=True)
args = vars(parser.parse_args())

input_table_path = args['inputtable']

# First, parse table
input_table = open(input_table_path, 'r')
input_table_rows = input_table.read().splitlines()
# Note - assume table is column 0 and product is column 4
table_products = {}
unique_markers = {}
for i,line in enumerate(input_table_rows):
	if i > 0:
		line_list = line.split('\t')
		table_products[line_list[0]] = float(line_list[4])
		unique_markers[line_list[0]] = int(line_list[1])

# Sort tables in descending order of product
sorted_tables = sorted(table_products, key=table_products.get, reverse=True)

# Select all tables that have the same product as the top item
selected_unique_markers = {}
selected_unique_markers[sorted_tables[0]] = unique_markers[sorted_tables[0]]
for i,table in enumerate(sorted_tables):
	if i > 0:
		if table_products[sorted_tables[i]] == table_products[sorted_tables[0]]:
			selected_unique_markers[sorted_tables[i]] = unique_markers[sorted_tables[i]]

# Sort all selected tables in descending order of number of unique markers
sorted_selected_tables = sorted(selected_unique_markers, key=selected_unique_markers.get, reverse=True)

print (sorted_selected_tables[0])
