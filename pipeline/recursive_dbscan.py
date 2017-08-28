#!/usr/bin/env python

import pandas as pd
from sklearn.cluster import DBSCAN
from scipy import stats
import sys
import copy
import numpy as np
import numbers
import math
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from tsne import bh_sne
from sklearn import decomposition
import os
import statistics
import argparse
import logging
import pdb


def run_BH_tSNE(table, do_pca=True):

	pca_dimensions = 50
	perplexity = 30.0

	logger.info("run_BH_tSNE: Running k-mer based binning...")
	# Note - currently doesn't handle cases where PCA dimensions and perplexity set too high

	# We make a submatrix, consisting of the contigs in the table
	k_mer_counts_submatrix = list()
	for i,row in table.iterrows():
		k_mer_counts_submatrix.append(k_mer_counts[i])

	normalized_k_mer_submatrix = normalizeKmers(k_mer_counts_submatrix)

	# PCA

	if (len(normalized_k_mer_submatrix[0]) > pca_dimensions) and (do_pca == True):
		logger.info('run_BH_tSNE: Principal component analysis')
		pca = decomposition.PCA(n_components=pca_dimensions)
		pca_matrix = pca.fit_transform(normalized_k_mer_submatrix)
	else:
		logger.info('run_BH_tSNE: Principle component analysis step skipped')

	# BH-tSNE
	logger.info('run_BH_tSNE: BH-tSNE')

	# Adjust perplexity according to the number of data points
	# Took logic from tsne source code
	if (len(normalized_k_mer_submatrix) - 1) < (3 * perplexity)  :
		perplexity = (float(len(normalized_k_mer_submatrix) - 1) / 3) - 1

	logger.info(str(len(normalized_k_mer_submatrix)) + ' data points')
	logger.info(str(len(normalized_k_mer_submatrix[0])) + ' dimensions')

	if (len(normalized_k_mer_submatrix[0]) > pca_dimensions) and (do_pca == True):
		X = np.array(pca_matrix)
	else:
		X = np.array(normalized_k_mer_submatrix)
	bh_tsne_matrix = bh_sne(X, d=2, perplexity=perplexity, theta=0.5)

	# We will add bh_tsne_x and bh_tsne_y columns to the contig table

	bh_tsne_x = list()
	bh_tsne_y = list()

	for i in range(0, len(bh_tsne_matrix)):
		bh_tsne_x.append(bh_tsne_matrix[i][0])
		bh_tsne_y.append(bh_tsne_matrix[i][1])

	table['bh_tsne_x'] = pd.Series(bh_tsne_x, index = table.index)
	table['bh_tsne_y'] = pd.Series(bh_tsne_y, index = table.index)

def runDBSCANs(table, dimensions):
	# Carry out DBSCAN, starting at eps=0.3 and continuing until there is just one group
	current_eps = 0.3
	db_tables = {} # Will be keyed by eps
	number_of_clusters = float('inf')
	current_step = 0.1
	number_of_tables = {}
	while(number_of_clusters > 1):
		logger.info('EPS: ' + str(current_eps))
		dbscan_output_pd = dbscan_simple(table, current_eps, dimensions)
		new_pd_copy = copy.deepcopy(dbscan_output_pd)
		db_tables[current_eps] = new_pd_copy
		current_eps = current_eps + current_step

		# Count the number of clusters
		number_of_clusters = countClusters(new_pd_copy)
		if number_of_clusters in number_of_tables:
			number_of_tables[number_of_clusters] += 1
		else:
			number_of_tables[number_of_clusters] = 1

		# We speed this up if we are getting a lot of tables with the same number of clusters
		if number_of_tables[number_of_clusters] > 10:
			current_step = current_step * 10

		logger.info('number of clusters: ' + str(number_of_clusters))

	return db_tables

def dbscan_simple(table, eps, dimensions):
	table_copy = copy.deepcopy(table)
	table_size = len(table.index)
	#logger.debug('dbscan_simple, eps: ' + str(eps) + ', table_size: ' + str(table_size))
	# Delete db_cluster column
	if 'db_cluster' in table_copy:
		table_copy.drop('db_cluster')

	# Make a matrix
	if dimensions == 2:
		X = table_copy.as_matrix(columns=['bh_tsne_x', 'bh_tsne_y'])
	elif dimensions == 3:
		X = table_copy.as_matrix(columns=['bh_tsne_x', 'bh_tsne_y', 'cov'])
	db = DBSCAN(eps=eps, min_samples=1).fit(X)

	table_copy['db_cluster'] = db.labels_

	return table_copy

def countClusters(pandas_table):
	clusters = {}
	for i, row in pandas_table.iterrows():
		cluster = row['db_cluster']
		if cluster not in clusters:
			clusters[cluster] = 1
		else:
			clusters[cluster] += 1
	number_of_clusters = len(list(clusters.keys()))
	return number_of_clusters

def getClusterInfo(pandas_table, hmm_dictionary, life_domain):
	marker_totals = {}
	for i, row in pandas_table.iterrows():
		contig = row['contig']
		cluster = row['db_cluster']
		if cluster not in marker_totals:
			marker_totals[cluster] = {}

		if contig not in hmm_dictionary:
			continue

		for pfam in hmm_dictionary[contig]:
			if pfam in marker_totals[cluster]:
				marker_totals[cluster][pfam] += hmm_dictionary[contig][pfam]
			else:
				marker_totals[cluster][pfam] = hmm_dictionary[contig][pfam]

	expected_number = 139
	if life_domain == 'archaea':
		expected_number = 162

	cluster_details = {} # Will hold completeness, purity

	for cluster in marker_totals:
		total_markers = len(marker_totals[cluster])
		total_marker_count = 0
		total_unique = 0
		for marker in marker_totals[cluster]:
			total_marker_count += marker_totals[cluster][marker]
			if marker_totals[cluster][marker] == 1:
				total_unique += 1
		completeness = (float(total_markers) / expected_number) * 100
		# Protect from divide by zero
		if total_markers == 0:
			purity = 0
		else:
			purity = (float(total_unique) / total_markers) * 100

		cluster_details[cluster] = { 'completeness': completeness, 'purity': purity }

	return cluster_details

def assessDBSCAN(table_dictionary, hmm_dictionary, domain, completeness_cutoff, purity_cutoff):
	# Assess clusters of each DBSCAN table
	cluster_info = {} # Dictionary that is keyed by eps, will hold details of each cluster in each table
	for eps in table_dictionary:
		current_table = table_dictionary[eps]
		current_cluster_info = getClusterInfo(current_table, hmm_dictionary, domain)
		cluster_info[eps] = current_cluster_info

	number_complete_and_pure_clusters = {}
	number_pure_clusters = {}
	median_completeness = {}
	mean_completeness = {}
	completeness_product = {}
	for eps in cluster_info:
		complete_clusters = 0
		pure_clusters = 0
		completenessList = []
		for cluster in cluster_info[eps]:
			completeness = cluster_info[eps][cluster]['completeness']
			purity = cluster_info[eps][cluster]['purity']
			if completeness > completeness_cutoff and purity > purity_cutoff:
				complete_clusters += 1
				completenessList.append(completeness)
			if purity > purity_cutoff:
				pure_clusters += 1
		number_complete_and_pure_clusters[eps] = complete_clusters
		number_pure_clusters[eps] = pure_clusters
		# Protect against warning if list is empty
		if completenessList:
			median_completeness[eps] = np.median(completenessList)
			mean_completeness[eps] = np.mean(completenessList)
		else:
			median_completeness[eps] = 0
			mean_completeness[eps] = 0

		completeness_product[eps] = number_complete_and_pure_clusters[eps] * mean_completeness[eps]

	# Get eps value with highest number of complete clusters
	#sorted_eps_values = sorted(number_complete_and_pure_clusters, key=number_complete_and_pure_clusters.__getitem__, reverse=True)
	sorted_eps_values = sorted(median_completeness, key=median_completeness.__getitem__, reverse=True)
	#sorted_eps_values = sorted(mean_completeness, key=mean_completeness.__getitem__, reverse=True)
	best_eps_value = sorted_eps_values[0]

	# For impure clusters, output BH_tSNE table
	# First, find pure clusters
	best_db_table = table_dictionary[best_eps_value]

	complete_and_pure_clusters = {}
	other_clusters = {}
	for cluster in cluster_info[best_eps_value]:
		completeness = cluster_info[best_eps_value][cluster]['completeness']
		purity = cluster_info[best_eps_value][cluster]['purity']

		# Explicitly remove the noise cluster (-1)
		# Note: in the R DBSCAN implementation, the noise cluster is 0, in the sklearn implementation, the noise cluster is -1
		if cluster == -1:
			other_clusters[cluster] = 1
			continue

		if completeness > completeness_cutoff and purity > purity_cutoff:
			complete_and_pure_clusters[cluster] = 1
		else:
			other_clusters[cluster] = 1

	# Subset the data frame
	subset_other_db_table = copy.deepcopy(best_db_table)

	for cluster in complete_and_pure_clusters:
		subset_other_db_table = subset_other_db_table[subset_other_db_table['db_cluster'] != cluster]

	# Output a table with the classifications, mainly for debug/investigation purposes
	output_db_table = copy.deepcopy(best_db_table)
	for index, row in output_db_table.iterrows():
		if row['db_cluster'] in other_clusters:
			output_db_table.set_value(index, 'db_cluster', -1)

	# We now drop the db_cluster column
	subset_other_db_table = subset_other_db_table.drop('db_cluster', 1)
	unclustered = subset_other_db_table

	# We now make a data structure containing cluster information for complete clusters only
	output_cluster_info = {}
	output_contig_cluster = {}
	for cluster in complete_and_pure_clusters:
		output_cluster_info[cluster] = {'completeness': cluster_info[best_eps_value][cluster]['completeness'], 'purity': cluster_info[best_eps_value][cluster]['purity']}

	# Now we grab contig names from the best db table
	for i, row in best_db_table.iterrows():
		contig = row['contig']
		cluster = row['db_cluster']
		if cluster in complete_and_pure_clusters:
			output_contig_cluster[contig] = cluster

	return output_cluster_info, output_contig_cluster, unclustered

def revcomp( string ):
	trans_dict = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' }
	complement_list = list()

	for i in range(0, len(string) ):
		if string[i] in trans_dict:
			complement_list.append(trans_dict[string[i]])
		else:
			return -1

	return ''.join(reversed(complement_list))

def normalizeKmers(count_matrix): # list of lists, not a np matrix
	# We now remove all the k-mers where all counts are '1'
	logger.info('Trimming k-mers')
	columns_to_delete = dict()
	for i in range(0, len(unique_k_mers.keys())):
		non_zero_counts = 0
		for j in range(0, len(count_matrix)):
			if count_matrix[j][i] > 1:
				non_zero_counts += 1

		if non_zero_counts == 0:
			columns_to_delete[i] = 1


	filtered_contig_k_mer_counts = list()
	for i in range(0, len(count_matrix)):
		new_row = list()
		for j in range(0, len(unique_k_mers.keys())):
			if j not in columns_to_delete:
				new_row.append(count_matrix[i][j])
		filtered_contig_k_mer_counts.append(new_row)

	# 3. Normalization
	logger.info('Normalizing counts')

	k_mer_frequency_matrix = list()

	for count_list in filtered_contig_k_mer_counts:
		total_count = 0
		for count in count_list:
			total_count += count

		normalized_list = list()
		for count in count_list:
			normalized_list.append(float(count)/total_count)

		# Now we calculate the Centered log-ratio (CLR) transformation
		# See Aitchison, J. The Statistical Analysis of Compositional Data (1986) and
		# Pawlowsky-Glahn, Egozcue, Tolosana-Delgado. Lecture Notes on Compositional Data Analysis (2011)
		geometric_mean = stats.mstats.gmean(normalized_list)
		clr_list = list()

		for item in normalized_list:
			intermediate_value = item / geometric_mean
			clr_list.append(math.log(intermediate_value))

		k_mer_frequency_matrix.append(clr_list)

	return k_mer_frequency_matrix

parser = argparse.ArgumentParser(description="Perform initial clustering via BH-tSNE and DBSCAN.")
parser.add_argument('-t','--input_table', help='Master contig table. Optionally can contain taxonomy data', required=True)
parser.add_argument('-a','--assembly_fasta', help='Assembly fasta', required=True)
parser.add_argument('-o','--output_table', help='Path to output table', required=True)
parser.add_argument('-k','--kingdom', help='Kingdom to consider (archaea|bacteria)', choices=['bacteria','archaea'], default = 'bacteria')

args = vars(parser.parse_args())

input_table_path = args['input_table']
input_fasta_path = args['assembly_fasta']
output_table_path = args['output_table']
domain = args['kingdom']

#logger
logger = logging.getLogger('recursive_dbscan.py')
hdlr = logging.FileHandler('recursive_dbscan.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
console = logging.StreamHandler()
console.setFormatter(formatter)
logger.setLevel(logging.DEBUG)
logger.addHandler(console)

master_table = pd.read_table(input_table_path)
master_table['cluster'] = 'unclustered'
master_table['bh_tsne_x'] = 0
master_table['bh_tsne_y'] = 0

# Load fasta sequences
assembly_seqs = {}
for seq_record in SeqIO.parse(input_fasta_path, 'fasta'):
	assembly_seqs[seq_record.id] = seq_record

# Count K-mer frequencies
k_mer_size = 5
matrix_file = 'k-mer_matrix'
k_mer_dict = dict() # Holds lists of k-mer counts, keyed by contig name

count = 0
unique_k_mers = dict()

DNA_letters = ['A', 'T', 'C', 'G']
all_k_mers = list(DNA_letters)
for i in range(1, k_mer_size):
	new_list = list()
	for current_seq in all_k_mers:
		for char in DNA_letters:
			new_list.append(current_seq + char)
	all_k_mers = new_list

# Now we trim k-mers and put them in the dictionary
for k_mer in all_k_mers:
	k_mer_reverse = revcomp(k_mer)
	if (k_mer not in unique_k_mers) and (k_mer_reverse not in unique_k_mers):
		unique_k_mers[k_mer] = count
		count += 1

if os.path.isfile(matrix_file):
	logger.info("K-mer matrix already exists!")
	logger.info("Continuing to next step...")

	# Now we load the k-mer matrix
	with open(matrix_file) as matrix:
		for i,line in enumerate(matrix):
			if i > 0:
				line_list = line.rstrip().split('\t')
				contig = line_list.pop(0)
				line_list = [ int(x) for x in line_list ]
				k_mer_dict[contig] = line_list
else:
	# Count k-mers
	# First we make a dictionary of all the possible k-mers (discounting revcomps)
	# Under each key is an index to be used in the subsequent lists
	# The order of the indices depends on the order k-mers were encountered while making the dictionary
	logger.info('Counting k-mers')

	for contig_name in assembly_seqs:
		contig_seq = str(assembly_seqs[contig_name].seq)
		# Initialize the list for the contig to be all 1s - as we can't have zero values in there for CLR later
		list_size = len(unique_k_mers.keys())
		current_contig_k_mer_counts = [1 for k in range(0, list_size)]

		for j in range(0, (len(contig_seq) - k_mer_size)):
			k_mer = contig_seq[j:j + k_mer_size]
			k_mer_reverse = revcomp(k_mer)

			# Find appropriate index
			# Note - this part naturally ignores any k_mers with weird characters
			if k_mer in unique_k_mers:
				index = unique_k_mers[k_mer]
				current_contig_k_mer_counts[index] += 1
			elif k_mer_reverse in unique_k_mers:
				index = unique_k_mers[k_mer_reverse]
				current_contig_k_mer_counts[index] += 1

		k_mer_dict[contig_name] = current_contig_k_mer_counts

	# Write the file in case we have to do this again
	matrix = open(matrix_file, 'w')

	# The first line consists of the k-mer headings (left corner is blank because the contig names are listed under it)
	header_line_list = [''] + sorted(unique_k_mers, key=unique_k_mers.__getitem__)
	header_line = '\t'.join(header_line_list)
	matrix.write(header_line + '\n')

	for contig in k_mer_dict:
		line_list = k_mer_dict[contig]
		output_list = [contig] + line_list
		for i, item in enumerate(output_list):
			output_list[i] = str(item)
		output_line = '\t'.join(output_list)
		matrix.write(output_line + '\n')
	matrix.close()

### Collate training data for ML steps later
# We now set up global data structures to be used in supervised machine learning
contig_list = master_table['contig'].tolist()
coverage_list = master_table['cov'].tolist()
taxonomy_matrix = list()

# Make normalized k-mer matrix
k_mer_counts = list()
for contig in contig_list:
	k_mer_counts.append(k_mer_dict[contig])

normalized_k_mer_matrix = normalizeKmers(k_mer_counts)




run_BH_tSNE(master_table)

contig_markers = {}
for i, row in master_table.iterrows():
    contig = row['contig']
    pfamString = row['single_copy_PFAMs']
    if pfamString == 'NA' or (isinstance(pfamString, numbers.Number) and math.isnan(pfamString)):
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



completeness_cutoff = 20
purity_cutoff = 90
round_counter = 0
global_cluster_info = {}
local_current_table = copy.deepcopy(master_table)

data_size = len(master_table.index)

if 'kingdom' in master_table.columns:
	has_taxonomy_info = True
else:
	has_taxonomy_info = False

if has_taxonomy_info and data_size > 50:
	for dimensions in [2, 3]:
		taxonomic_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
		logger.info('Further splitting according to taxonomic classifications')
		for taxonomic_level in taxonomic_levels:
			logger.info('Taxonomic level: ' + taxonomic_level)
			unclustered_table = pd.DataFrame()

			# Make subtables for each type of classification at the current level
			classification_dict = dict()
			for i,row in local_current_table.iterrows():
				classification_dict[row[taxonomic_level]] = 1

			# Skip iteration if the current taxonomic level is empty
			if not classification_dict:
				continue

			for classification in classification_dict.keys():
				logger.info('Examining ' + classification)
				# Get subset table
				subset_table = local_current_table.loc[local_current_table[taxonomic_level] == classification] # COPY of local_current_table

				while True:
					if len(subset_table.index) < 1:
						break
					round_counter += 1
					logger.info('Running DBSCAN round ' + str(round_counter))

					db_tables = runDBSCANs(subset_table, dimensions)
					cluster_information, contig_cluster_dictionary, local_unclustered_table = assessDBSCAN(db_tables, contig_markers, domain, completeness_cutoff, purity_cutoff)

					subset_table = local_unclustered_table

					if not cluster_information:
						break

					# Populate the global data structures
					for	cluster in cluster_information:
						new_cluster_name = 'DBSCAN' + '_round' + str(round_counter) + '_' + str(cluster)
						global_cluster_info[new_cluster_name] = cluster_information[cluster]

					for contig in contig_cluster_dictionary:
						new_cluster_name = 'DBSCAN'+ '_round' + str(round_counter) + '_' + str(contig_cluster_dictionary[contig])
						table_indices = local_current_table[local_current_table['contig'] == contig].index.tolist()
						master_table.set_value(table_indices[0], 'cluster', new_cluster_name)

				# Add unclustered_table to combined unclustered dataframe
				unclustered_table = unclustered_table.append(local_unclustered_table)

			#current_table = copy.deepcopy(unclustered_table)
			local_current_table = unclustered_table
else:
	for dimensions in [2, 3]:
		while True:
		    round_counter += 1
		    logger.info('Running DBSCAN round ' + str(round_counter))

		    db_tables = runDBSCANs(local_current_table, dimensions)
		    cluster_information, contig_cluster_dictionary, unclustered_table = assessDBSCAN(db_tables, contig_markers, domain, completeness_cutoff, purity_cutoff)

		    if not cluster_information:
		        break

		    # Populate the global data structures
		    for	cluster in cluster_information:
		        new_cluster_name = 'DBSCAN' + '_round' + str(round_counter) + '_' + str(cluster)
		        global_cluster_info[new_cluster_name] = cluster_information[cluster]

		    for contig in contig_cluster_dictionary:
		        new_cluster_name = 'DBSCAN' + '_round' + str(round_counter) + '_' + str(contig_cluster_dictionary[contig])
		        table_indices = master_table[master_table['contig'] == contig].index.tolist()
		        master_table.set_value(table_indices[0], 'cluster', new_cluster_name)

		    local_current_table = unclustered_table

# Output table
master_table.to_csv(path_or_buf=output_table_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)
