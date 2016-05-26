#!/usr/bin/env python

import sys
import argparse 
import pdb
import numpy

parser = argparse.ArgumentParser(description='Script to assess the metagenome deconvolution using single-copy marker genes')
parser.add_argument('-s','--hmmtable', help='HMM table created by make_marker_table.py', required=True)
parser.add_argument('-d', '--dbscantable', help='Table(s) containing DBSCAN clusters', nargs='+', required=True)
parser.add_argument('-o', '--output', help='Output table path', required=True)
parser.add_argument('-c', '--column', help='bin column name', default = 'db.cluster')
parser.add_argument('-k', '--kingdom', help='Kingdom under consideration (bacteria|archaea)', default = 'bacteria')
args = vars(parser.parse_args())

hmm_contig_table_path = args['hmmtable']
dbscan_table_paths = args['dbscantable']
output_table_path = args['output']
dbscan_column_heading = args['column']
kingdom = args['kingdom']

# First parse hmm_contig_table_path

hmm_contig_table_rows = None
with open(hmm_contig_table_path) as hmm_contig_table:
	hmm_contig_table_rows = hmm_contig_table.read().splitlines()

contig_markers = {}
for i,line in enumerate(hmm_contig_table_rows):
	if i > 0:
		line_list = line.split('\t')
		pfam_list = line_list[1].split(',')
		contig = line_list[0]
		contig_markers[contig] = {}
		for pfam in pfam_list:
			if pfam in contig_markers[contig]:
				contig_markers[contig][pfam] += 1
			else:
				contig_markers[contig][pfam] = 1

# Now we go through the dbscan tables
#table_completeness_averages = {}
#table_contamination_averages = {}
table_binned_unique_marker_counter = {} # Keeps a total of unique markers in each bin, as long as the bin contains more than 20% of the total
# In Bacteria, the total is 139, in Archaea, the total is 162
table_numbers_of_clusters = {}


for dbscan_table_path in dbscan_table_paths:
	print ('Considering ' + dbscan_table_path)

	bin_markers = {}
	dbscan_table_rows = None
	with open(dbscan_table_path) as dbscan_table:
		dbscan_table_rows = dbscan_table.read().splitlines()

	bin_column_index = None
	number_found = 0
	first_line_list = dbscan_table_rows[0].split('\t')
	for i,value in enumerate(first_line_list):
		if value == dbscan_column_heading:
			bin_column_index = i
			number_found += 1
	if number_found > 1:
		print 'Error, bin table has more than one column headed ' + dbscan_column_heading
		sys.exit(2)
	if not bin_column_index:
		print 'Error, could not find column ' + dbscan_column_heading + ' in bin table ' + dbscan_table_path
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

	for i, line in enumerate(dbscan_table_rows):
		if i > 0:
			line_list = line.split('\t')
			contig = line_list[contig_column_index]
			cluster = line_list[bin_column_index]

			if cluster not in bin_markers:
				bin_markers[cluster] = {}

			if contig in contig_markers:
				for pfam in contig_markers[contig]:
					if pfam in bin_markers[cluster]:
						bin_markers[cluster][pfam] += contig_markers[contig][pfam]
					else:
						bin_markers[cluster][pfam] = contig_markers[contig][pfam]

	# Now work out the number duplicated for each cluster
	#duplicated = [] # Just a list of counts for duplicated markers in each cluster
	table_unique_markers_counter = 0
	number_of_clusters_over_threshold = 0
	for cluster in bin_markers:
		not_duplicated_counter = 0
		for pfam in bin_markers[cluster]:
			if bin_markers[cluster][pfam] == 1:
				not_duplicated_counter += 1
		#duplicated.append(duplicated_counter)
		threshold = None
		if kingdom == 'bacteria':
			threshold = 28
		elif kingdom == 'archaea':
			threshold = 32
		if not_duplicated_counter > threshold:
			table_unique_markers_counter += not_duplicated_counter
			number_of_clusters_over_threshold += 1

	table_binned_unique_marker_counter[dbscan_table_path] = table_unique_markers_counter
	table_numbers_of_clusters[dbscan_table_path] = number_of_clusters_over_threshold

	#average_duplicated = sum(duplicated) / len(duplicated)
	#table_contamination_averages[dbscan_table_path] = average_duplicated

	# Now work out the number of markers each cluster has
	#number_of_markers = [] # Just a list of counts of unique markers in each cluster
	#for cluster in bin_markers:
	#	number_of_markers.append(len(bin_markers[cluster]))

	#average_markers = sum(number_of_markers) / len(number_of_markers)
	#table_completeness_averages[dbscan_table_path] = average_markers


# Print output table
output_table = open(output_table_path, 'w')
#output_table.write('table\tav_number_of_markers\tav_duplicated_markers\n')
output_table.write('table\tnumber_binned_unique_markers\tnumber_of_clusters_over_threshold\n')

#for table_path in table_completeness_averages:
	#completeness = table_completeness_averages[table_path]
	#duplicated = table_contamination_averages[table_path]
	#output_table.write(table_path + '\t' + str(completeness) + '\t' + str(duplicated) + '\n')
for table_path in table_binned_unique_marker_counter:
	number_unique_markers = table_binned_unique_marker_counter[table_path]
	number_of_clusters = table_numbers_of_clusters[table_path]
	output_table.write(table_path + '\t' + str(number_unique_markers) + '\t' + str(number_of_clusters) + '\n')