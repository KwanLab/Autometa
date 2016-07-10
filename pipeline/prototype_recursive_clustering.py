#!/usr/bin/env python

# Program to carry out recursive clustering using vizbin coordinates using machine learning techniques in R

import rpy2.robjects as robjects # Bridge to R code
from rpy2.robjects.packages import importr
import sys
import os
import pandas as pd
import csv
from rpy2.robjects import pandas2ri
import numpy
import math
from Bio import SeqIO
import argparse
import pprint
import pdb

pandas2ri.activate()

# Load R libaries
rbase = importr('base')
dbscan = importr('dbscan')
clv = importr('clv')
LICORS = importr('LICORS')

# Define R functions
robjects.r('''
	dbscan_run <- function(path) {
		input_data <- read.table(path, header=TRUE)
		# Funny things happen if there is already a 'db_cluster' column, let's delete it!
		if ("db_cluster" %in% colnames(input_data))
		{
			input_data$db_cluster <- NULL
		}

		d <- data.frame(input_data$vizbin_x, input_data$vizbin_y)

		# Cycle through eps values 0.3 to 5.0
		eps_values = seq(from = 0.3, to = 5.0, by = 0.1)
		db_values = vector(mode = "numeric", length = length(eps_values))

		for (i in seq(length.out=length(eps_values))) {
			db <- dbscan(d, eps=eps_values[i], minPts=3)
			full_table <- data.frame(input_data, db_cluster=db$cluster)
			full_table_no_noise <- subset(full_table, db_cluster != '0')
			# Count number of unique groups
			freqTable <- as.data.frame(table(full_table_no_noise$db_cluster))
			if (nrow(freqTable) == 1) {
				db_values[i] = Inf
				next
			}
			d2 <- data.frame(full_table_no_noise$vizbin_x, full_table_no_noise$vizbin_y)
			db2 <- as.vector(full_table_no_noise$db_cluster)
			index.list <- cls.scatt.data(d2, db2, dist="euclidean")
			davies.bouldin <- clv.Davies.Bouldin( index.list, "centroid", "centroid")
			db_values[i] <- davies.bouldin
		}

		# Decide which eps value is the best
		db_data_frame <- data.frame(dbscan_eps = eps_values, davies_bouldin_index = db_values)
		print(db_data_frame)
		db_sorted_data_frame <- db_data_frame[order(db_data_frame$davies_bouldin_index),]
		smallest_db_value <- db_sorted_data_frame[1:1, 2:2]
		best_eps_values <- c(db_sorted_data_frame[1:1, 1:1])

		for (i in seq(from=2, to=length(db_sorted_data_frame))) {
			if (db_sorted_data_frame[i:i, 2:2] == smallest_db_value) {
				best_eps_values <- c(best_eps_values, db_sorted_data_frame[i:i, 1:1])
			}
		}

		# Find median, to nearest 0.1
		best_eps_value <- round(median(best_eps_values), digits=1)
		print(paste('Best EPS value: ', best_eps_value))

		# Redo dbscan for the best eps value
		db <- dbscan(d, eps=best_eps_value, minPts=3)
		output_table <- data.frame(input_data, db_cluster = db$cluster )
		# print ( output_table )
		return( output_table )
	}

	kmeanspp_run <- function(input_data_frame, k_upper_bound) {
		if ("kmeanspp_cluster" %in% colnames(input_data_frame))
		{
			input_data_frame$kmeanspp_cluster <- NULL
		}

		d <- data.frame(input_data_frame$vizbin_x, input_data_frame$vizbin_y)
		k_upper_bound <- ceiling(k_upper_bound)

		if (k_upper_bound == 2) {
			# There is just one possible value of k, so let's just run one kmeanspp round
			cluster_result <- kmeanspp(d, k = k_upper_bound, start="random")
			full_table <- data.frame(input_data_frame, kmeanspp_cluster = cluster_result$cluster)
			return(full_table)
		}

		if (nrow(input_data_frame) == 2) {
			full_table <- data.frame(input_data_frame, kmeanspp_cluster = c(1, 2))
			return(full_table)
		}

		k_values = seq(from = 2, to = k_upper_bound)
		

		if ((nrow(input_data_frame) - 1) < k_upper_bound) {
			highest_i = nrow(input_data_frame) - 2
		} else {
			highest_i = length(k_values)
		}
		edited_k_values = vector(mode="numeric", length=length(seq(highest_i)))
		for ( i in seq(highest_i) ) {
			edited_k_values[i] = k_values[i]
		}
		db_values = vector(mode = "numeric", length = length(edited_k_values))

		for ( i in seq(highest_i) ) {
			current_k = edited_k_values[i]
			cluster_result <- kmeanspp(d, k = current_k, start = "random")
			full_table <- data.frame(input_data_frame, kmeanspp_cluster = cluster_result$cluster )
			cluster <- as.vector(full_table$kmeanspp_cluster)
			index.list <- cls.scatt.data(d, cluster, dist="euclidean")
			davies.bouldin <- clv.Davies.Bouldin(index.list, "centroid", "centroid")
			db_values[i] <- davies.bouldin
		}

		db_table <- data.frame(k_values = edited_k_values, davies_bouldin_index = db_values)
		db_sorted_data_frame <- db_table[order(db_table$davies_bouldin_index),]
		smallest_db_value <- db_sorted_data_frame[1:1, 2:2]
		best_k_values <- c(db_sorted_data_frame[1:1, 1:1])

		for (i in seq(from=2, to=length(db_sorted_data_frame))) {
			if (db_sorted_data_frame[i:i, 2:2] == smallest_db_value) {
				best_k_values <- c(best_k_values, db_sorted_data_frame[i:i, 1:1])
			}
		}

		# Find median, to nearest 1
		best_k_value <- floor(median(best_k_values))
		print(paste('Best K value: ', best_k_value))

		# Redo kmeanspp for best k value
		cluster_result <- kmeanspp(d, k = best_k_value, start = "random")
		full_table <- data.frame(input_data_frame, kmeanspp_cluster = cluster_result$cluster )
	
		return(full_table)
	}
	''')

def assess_cluster_completeness(pandas_table, marker_dictionary, life_domain):
	marker_totals = {}
	for i, row in pandas_table.iterrows():
		contig = row['contig']
		if contig in marker_dictionary:
			for marker in marker_dictionary[contig]:
				if marker in marker_totals:
					marker_totals[marker] += marker_dictionary[contig][marker]
				else:
					marker_totals[marker] = marker_dictionary[contig][marker]

	# Now work out the completeness and purity
	expected_number = 139
	if life_domain == 'archaea':
		expected_number = 162

	total_markers = len(marker_totals)
	total_unique = 0

	for marker in marker_totals:
		if marker_totals[marker] == 1:
			total_unique += 1

	completeness = (float(total_markers) / expected_number) * 100
	purity = (float(total_unique) / expected_number) * 100

	return(completeness, purity)

def upperBound(pandas_table, marker_dictionary):
	marker_totals = {}
	for i, row in pandas_table.iterrows():
		contig = row['contig']
		if contig in marker_dictionary:
			for marker in marker_dictionary[contig]:
				if marker in marker_totals:
					marker_totals[marker] += marker_dictionary[contig][marker]
				else:
					marker_totals[marker] = marker_dictionary[contig][marker]

	total_list = []
	for marker in marker_totals:
		total_list.append(marker_totals[marker])

	mean = numpy.mean(total_list)
	sd = numpy.std(total_list)

	upper_bound = mean + (3*sd)
	return upper_bound

def assessClusters(pandas_table, marker_dictionary, life_domain):
	# First get info on clusters
	clusters_present = {}
	for i, row in pandas_table.iterrows():
		cluster = row['kmeanspp_cluster']
		clusters_present[cluster] = 1

	cluster_info = {}
	cluster_completeness_info = {}
	cluster_purity_info = {}
	
	for cluster in clusters_present:
		subset_table = pandas_table.loc[(pandas_table.kmeanspp_cluster == cluster)]
		completeness, purity = assess_cluster_completeness(subset_table, marker_dictionary, life_domain)

		cluster_completeness_info[cluster] = completeness
		cluster_purity_info[cluster] = purity

		#print('cluster: ' + str(cluster) + ', completeness: ' + str(completeness) + ', purity: ' + str(purity))
		if completeness > 20 and purity > 90:
			cluster_info[cluster] = 'complete'
		elif completeness < 20:
			cluster_info[cluster] = 'failed'
		else:
			cluster_info[cluster] = 'impure'

	return cluster_info, cluster_completeness_info, cluster_purity_info

def getClusterInfo(pandas_table):
	cluster_contig_info = {} # Dictionary of dictionaries, keyed by cluster then by contig
	for i, row in pandas_table.iterrows():
		if row['db_cluster'] in cluster_contig_info:
			cluster_contig_info[row['db_cluster']][row['contig']] = { 'length': row['length'], 'gc': row['gc'], 'cov': row['cov']}
		else:
			cluster_contig_info[row['db_cluster']] = { row['contig']: { 'length': row['length'], 'gc': row['gc'], 'cov': row['cov']}}

	# Need to calculate weighted average of gc and cov, as well as N50
	cluster_contig_lengths = {} # Dictionary holding sorted lists of contig length (descending order)
	cluster_total_lengths = {} # Dictionary to hold total lengths

	for cluster in cluster_contig_info:
		cluster_contig_lengths[cluster] = []
		cluster_total_lengths[cluster] = 0
		for contig in cluster_contig_info[cluster]:
			length = cluster_contig_info[cluster][contig]['length']
			cluster_total_lengths[cluster] += length
			if not cluster_contig_lengths[cluster]: # List is empty
				cluster_contig_lengths[cluster] = [length]
			else:
				# Work out where in the list to put the new length
				insertion_index = None
				for i in range(len(cluster_contig_lengths[cluster])):
					if cluster_contig_lengths[cluster][i] < length:
						insertion_index = i
						break
				if insertion_index is None:
					cluster_contig_lengths[cluster].append(length)
				else:
					#print ('insertion_index: ' + str(i))
					cluster_contig_lengths[cluster].insert(insertion_index, length)

	cluster_n50s = {}
	for cluster in cluster_contig_lengths:
		target_length = float(cluster_total_lengths[cluster])/2
		running_total = 0
		n50 = None
		for current_length in cluster_contig_lengths[cluster]:
			running_total += current_length
			if running_total >= target_length:
				n50 = current_length
				break
		cluster_n50s[cluster] = n50

	# Need to calculate length fractions for weighted averages
	for cluster in cluster_contig_info:
		for contig in cluster_contig_info[cluster]:
			length = cluster_contig_info[cluster][contig]['length']
			length_fraction = float(length) / cluster_total_lengths[cluster]
			cluster_contig_info[cluster][contig]['length_fraction'] = length_fraction

	cluster_gc_weighted_av = {}
	cluster_cov_weighted_av = {}

	for cluster in cluster_contig_info:
		cluster_gc_weighted_av[cluster] = 0
		cluster_cov_weighted_av[cluster] = 0
		for contig in cluster_contig_info[cluster]:
			gc_addition = float(cluster_contig_info[cluster][contig]['gc']) * cluster_contig_info[cluster][contig]['length_fraction']
			cluster_gc_weighted_av[cluster] += gc_addition
			cov_addition = float(cluster_contig_info[cluster][contig]['cov']) * cluster_contig_info[cluster][contig]['length_fraction']
			cluster_cov_weighted_av[cluster] += cov_addition

	# Make the output data structure
	output_dictionary = {}
	for cluster in cluster_contig_info:
		output_dictionary[cluster] = { 'size': cluster_total_lengths[cluster], 'longest_contig': cluster_contig_lengths[cluster][0], 'n50': cluster_n50s[cluster], 'number_contigs': len(cluster_contig_lengths[cluster]), 'cov': cluster_cov_weighted_av[cluster], 'gc_percent': cluster_gc_weighted_av[cluster] }

	return output_dictionary

dbscan_run = robjects.r['dbscan_run']
kmeanspp_run = robjects.r['kmeanspp_run']

pp = pprint.PrettyPrinter(indent=4)

parser = argparse.ArgumentParser(description="Prototype script to automatically cluster using DBSCAN, Kmeans++ and Davies-Bouldin Index")
parser.add_argument('-m','--marker_tab', help='Output of make_marker_table.py', required=True)
parser.add_argument('-v','--vizbin_tab', help='Table containing vizbin coordinates', required=True)
parser.add_argument('-d','--domain', help='Microbial domain (bacteria|archaea)', default='bacteria')
parser.add_argument('-f','--fasta', help='Assembly FASTA file') # Optional, if present, will output cluster and nonclustered contigs
parser.add_argument('-o','--outdir', help='Path of directory for output', required=True)
args = vars(parser.parse_args())

hmm_table_path = args['marker_tab']
vizbin_table_path = args['vizbin_tab']
domain = args['domain']
fasta_path = args['fasta']
outdir = os.path.abspath(args['outdir'])

# 1. Parse hmm table
contig_markers = {}
hmm_table = open(hmm_table_path, 'r')
hmm_table_lines = hmm_table.read().splitlines()

for i, line in enumerate(hmm_table_lines):
	if i > 0:
		lineList = line.split('\t')
		contig = lineList[0]
		pfamString = lineList[1]
		if pfamString == 'NA':
			continue
		pfamList = pfamString.split(',')
		# Note: we assume here that each contig only occurs on one line in the table
		for pfam in pfamList:
			if contig in contig_markers:
				if pfam in contig_markers[contig]:
					contig_markers[contig][pfam] += 1
				else:
					contig_markers[contig][pfam] = 1
			else:
				contig_markers[contig] = { pfam: 1 }

# 2. Carry out initial dbscan clustering
# 2a. parse vizbin table
vizbin_table = open(vizbin_table_path, 'r')
vizbin_table_lines = vizbin_table.read().splitlines()
headerLine = vizbin_table_lines[0]
headerLineList = headerLine.split('\t')

# Do dbscan
abs_vizbin_path = os.path.abspath(vizbin_table_path)
# Run dbscan_run R function
dbscan_output = dbscan_run(abs_vizbin_path)

# Now we extract cluster information from the dataframe
dbscan_output_pd = pandas2ri.ri2py(dbscan_output)
dbscan_output_pd['db_cluster'].apply(str)
clusters = {}
cluster_completeness = {}
cluster_purity = {}

for i, row in dbscan_output_pd.iterrows():
	cluster = row['db_cluster']
	if cluster not in clusters:
		if cluster == '0':
			clusters[cluster] = 'noise'
		else:
			clusters[cluster] = 'awaiting_assessment'

# For each cluster, assess completeness
for cluster in clusters:
	if clusters[cluster] != 'noise':
		subset_table = dbscan_output_pd.loc[(dbscan_output_pd.db_cluster == cluster)]
		completeness, purity = assess_cluster_completeness(subset_table, contig_markers, domain)
		cluster_completeness[cluster] = completeness
		cluster_purity[cluster] = purity

		if completeness > 20 and purity > 90:
			clusters[cluster] = 'complete'
		elif completeness < 20:
			clusters[cluster] = 'failed'
		else:
			clusters[cluster] = 'impure'

# 3. Carry out subclustering (like recursive sort)

# Count number of impure clusters
impure_clusters = 0
for cluster in clusters:
	if clusters[cluster] == 'impure':
		impure_clusters += 1

while (impure_clusters):
	#print ('Impure clusters: ' + str(impure_clusters))
	clusters_to_consider = list(clusters.keys())
	for cluster in clusters_to_consider:
		if clusters[cluster] == 'impure':
			print ('Considering cluster: ' + str(cluster))
			#pdb.set_trace()
			subset_table = dbscan_output_pd.loc[(dbscan_output_pd.db_cluster == cluster)]
			#print(subset_table)
			# Determine upper bound
			upper_bound = upperBound(subset_table, contig_markers)
			if math.isnan(upper_bound) or upper_bound < 2:
				clusters[cluster] = 'failed'
				continue

			print ('upper bound: ' + str(upper_bound))
			# Convert to R dataframe
			subset_r = pandas2ri.py2ri(subset_table)
			pp.pprint(subset_table)
			new_clustering_r = kmeanspp_run(subset_r, upper_bound)
			new_clustering_pd = pandas2ri.ri2py(new_clustering_r)

			print ('subclustering:')
			#print (new_clustering_pd)
			new_cluster_info, new_cluster_completeness, new_cluster_purity = assessClusters(new_clustering_pd, contig_markers, domain)
			#pp.pprint(new_cluster_info)

			# If one impure group, then change to failed
			new_cluster_list = list(new_cluster_info.keys())
			if len(new_cluster_list) == 1:
				new_cluster_info[new_cluster_list[0]] = 'failed'

			# Delete old key in clusters
			clusters.pop(cluster, None)
			cluster_completeness.pop(cluster, None)
			cluster_purity.pop(cluster, None)

			# Add entries for new clusters
			for new_cluster in new_cluster_info:
				assigned_cluster_name = str(cluster) + '_' + str(new_cluster)
				clusters[assigned_cluster_name] = new_cluster_info[new_cluster]
				cluster_completeness[assigned_cluster_name] = new_cluster_completeness[new_cluster]
				cluster_purity[assigned_cluster_name] = new_cluster_purity[new_cluster]

			# Change entries in pandas table
			contig_clusters = {}
			for i, row in new_clustering_pd.iterrows():
				contig = row['contig']
				new_cluster = row['kmeanspp_cluster']
				assigned_cluster_name = str(cluster) + '_' + str(new_cluster)
				contig_clusters[contig] = assigned_cluster_name

			for i, row in dbscan_output_pd.iterrows():
				contig = row['contig']
				if contig in contig_clusters:
					dbscan_output_pd.ix[i, 'db_cluster'] = contig_clusters[contig]

	impure_clusters = 0
	# Count impure clusters again
	for cluster in clusters:
		if clusters[cluster] == 'impure':
			impure_clusters += 1

# Make summary table
cluster_info = getClusterInfo(dbscan_output_pd)

# print summary table
summary_table_path = outdir + '/summary_table'
summary_table = open(summary_table_path, 'w')
summary_table.write('cluster\tsize\tlongest_contig\tn50\tnumber_contigs\tcompleteness\tpurity\tcov\tgc_percent\tstatus\n')
for cluster in cluster_info:
	size = str(cluster_info[cluster]['size'])
	longest_contig = str(cluster_info[cluster]['longest_contig'])
	n50 = str(cluster_info[cluster]['n50'])
	number_contigs = str(cluster_info[cluster]['number_contigs'])
	completeness = str(cluster_completeness[cluster])
	purity = str(cluster_purity[cluster])
	cov = str(cluster_info[cluster]['cov'])
	gc = str(cluster_info[cluster]['gc_percent'])
	status = str(clusters[cluster])

	outputString = '\t'.join([str(cluster), size, longest_contig, n50, number_contigs, completeness, purity, cov, gc, status]) + '\n'

	summary_table.write(outputString)
summary_table.close()

# Output whole pandas table
pandas_output_path = outdir + '/full_clustering_table'
dbscan_output_pd.to_csv(path_or_buf=pandas_output_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)

# Now for each 'complete' cluster, we output a separate fasta file.  We collect all contigs from the 'failed' and 'noise' clusters and 
# output them to one fasta (unclustered.fasta)
if fasta_path:
	assembly_seqs = {}
	for seq_record in SeqIO.parse(fasta_path, 'fasta'):
		assembly_seqs[seq_record.id] = seq_record

	# Initialize output fasta lists
	output_fastas = {} # Keyed by cluster
	output_fastas['unclustered'] = []
	for cluster in clusters:
		if clusters[cluster] == 'complete':
			output_fastas[cluster] = []

	for i, row in dbscan_output_pd.iterrows():
		contig = row['contig']
		cluster = row['db_cluster']
		seq_record = assembly_seqs[contig]
		if clusters[cluster] == 'complete':
			output_fastas[cluster].append(seq_record)
		else:
			output_fastas['unclustered'].append(seq_record)

	# Now write fasta files
	for cluster in output_fastas:
		if cluster == 'unclustered':
			fasta_output_path = outdir + '/' + str(cluster) + '.fasta'
		else:
			fasta_output_path = outdir + '/cluster_' + str(cluster) + '.fasta'
		SeqIO.write(output_fastas[cluster], fasta_output_path, 'fasta')


