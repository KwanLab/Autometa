#!/usr/bin/env python

# Program to carry out secondary dbscan clustering on a re-run of BH_tSNE on unclustered contigs (using low "perplexity" setting)
# When you run BH_tSNE this way, Davies-Bouldin index seems to fail because the groups are quite diffuse, although DBSCAN seems to do
# fine with the correct eps value.  Here we judge DBSCAN results based on the number of pure clusters (sum of total purity)

import readline
#import rpy2.robjects as robjects # Bridge to R code
#from rpy2.robjects.packages import importr
#from rpy2.robjects import pandas2ri
import sys
import os
import pandas as pd
import csv
import argparse
import copy
import numpy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from sklearn import decomposition
from sklearn.cluster import DBSCAN
from scipy import stats
import math
from tsne import bh_sne
import numpy as np
import logging
import subprocess
import getpass
import time 
import multiprocessing
import pprint
import pdb

#logger
logger = logging.getLogger('recursive_dbscan.py')
hdlr = logging.FileHandler('recursive_dbscan.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.DEBUG)

pp = pprint.PrettyPrinter(indent=4)

def dbscan_simple(table, eps):
	table_copy = copy.deepcopy(table)
	# Delete db_cluster column
	if 'db_cluster' in table_copy:
		table_copy.drop('db_cluster')

	# Make a matrix
	X = table_copy.as_matrix(columns=['bh_tsne_x', 'bh_tsne_y'])
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

		cluster_multiple = float(total_marker_count) / total_markers
		cluster_details[cluster] = { 'completeness': completeness, 'purity': purity, 'multiple': cluster_multiple }

	return cluster_details

def getClusterSummaryInfo(pandas_table):
	cluster_contig_info = {} # Dictionary of dictionaries, keyed by cluster then by contig
	for i, row in pandas_table.iterrows():
		if row['cluster'] in cluster_contig_info:
			cluster_contig_info[row['cluster']][row['contig']] = { 'length': row['length'], 'gc': row['gc'], 'cov': row['cov']}
		else:
			cluster_contig_info[row['cluster']] = { row['contig']: { 'length': row['length'], 'gc': row['gc'], 'cov': row['cov']}}

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

def runDBSCANs(table):
	# Carry out DBSCAN, starting at eps=0.3 and continuing until there is just one group
	current_eps = 0.3
	db_tables = {} # Will be keyed by eps
	number_of_clusters = float('inf')
	while(number_of_clusters > 1):
		#dbscan_output_r = dbscan_simple(r_table, current_eps)
		#dbscan_output_pd = pandas2ri.ri2py(dbscan_output_r)
		dbscan_output_pd = dbscan_simple(table, current_eps)
		new_pd_copy = copy.deepcopy(dbscan_output_pd)
		db_tables[current_eps] = new_pd_copy
		current_eps = current_eps + 0.1

		# Count the number of clusters
		number_of_clusters = countClusters(new_pd_copy)

	return db_tables

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
			median_completeness[eps] = numpy.median(completenessList)
			mean_completeness[eps] = numpy.mean(completenessList)
		else:
			median_completeness[eps] = 0
			mean_completeness[eps] = 0

		completeness_product[eps] = number_complete_and_pure_clusters[eps] * mean_completeness[eps]

	# Get eps value with highest number of complete clusters
	#sorted_eps_values = sorted(number_complete_and_pure_clusters, key=number_complete_and_pure_clusters.__getitem__, reverse=True)
	#sorted_eps_values = sorted(median_completeness, key=median_completeness.__getitem__, reverse=True)
	sorted_eps_values = sorted(mean_completeness, key=mean_completeness.__getitem__, reverse=True)
	best_eps_value = sorted_eps_values[0]

	#pdb.set_trace()

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

	output_filename = 'BH_tSNE' + str(BH_tSNE_counter) + '_DBSCAN' + str(round_counter) + '_eps' + str(best_eps_value) + '.tab'
	output_db_table.to_csv(path_or_buf=output_filename, sep='\t', index=False, quoting=None)

	# We now make an r table from subset_other_db_table, with the db_cluster column stripped
	subset_other_db_table = subset_other_db_table.drop('db_cluster', 1)
	#unclustered_r = pandas2ri.py2ri(subset_other_db_table)
	unclustered = subset_other_db_table

	#subset_complete_db_table = best_db_table
	#for cluster in other_clusters:
	#	subset_complete_db_table = subset_complete_db_table[subset_complete_db_table['db_cluster'] != cluster]

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

def run_BH_tSNE(fasta, output_filename,contig_table_path, do_pca=True):
	k_mer_size = 5
	pca_dimensions = 50
	perplexity = 30.0

	if os.path.isfile(output_filename):
		print "BH_tSNE output already exists!"
		print "Continuing to next step..."
		logger.info('BH_tSNE output already exists!')
		logger.info('Continuing to next step...')
		return None
	else:
		print "run_BH_tSNE: Running k-mer based binning..."
		logger.info('run_BH_tSNE: Running k-mer based binning...')
		# Note - currently doesn't handle cases where PCA dimensions and perplexity set too high

		# 1. Load fasta
		sequences = list()
		sequence_names = dict()
		seq_counter = 0
		for seq_record in SeqIO.parse(fasta, 'fasta'):
			sequences.append(seq_record)
			seq_name = str(seq_record.id)
			sequence_names[seq_name] = seq_counter
			seq_counter += 1

		# 2. Count k-mers
		# First we make a dictionary of all the possible k-mers (discounting revcomps)
		# Under each key is an index to be used in the subsequent lists
		# The order of the indices depends on the order k-mers were encountered while making the dictionary
		logger.info('run_BH_tSNE: Counting k-mers')
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

		contig_k_mer_counts = list()

		for i in range(0, len(sequences)):
			contig_name = str(sequences[i].id)
			contig_seq = str(sequences[i].seq)
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

			contig_k_mer_counts.append(current_contig_k_mer_counts)

		# We now remove all the k-mers where all counts are '1'
		logger.info('run_BH_tSNE: Trimming k-mers')
		columns_to_delete = dict()
		for i in range(0, len(unique_k_mers.keys())):
			non_zero_counts = 0
			for j in range(0, len(contig_k_mer_counts)):
				if contig_k_mer_counts[j][i] > 1:
					non_zero_counts += 1

			if non_zero_counts == 0:
				columns_to_delete[i] = 1


		filtered_contig_k_mer_counts = list()
		for i in range(0, len(contig_k_mer_counts)):
			new_row = list()
			for j in range(0, len(unique_k_mers.keys())):
				if j not in columns_to_delete:
					new_row.append(contig_k_mer_counts[i][j])
			filtered_contig_k_mer_counts.append(new_row)
		filtered_contig_k_mer_counts = contig_k_mer_counts

		# 3. Normalization
		logger.info('run_BH_tSNE: Normalizing counts')

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

		# 4. PCA

		if (len(k_mer_frequency_matrix[0]) > pca_dimensions) and (do_pca == True):
			logger.info('run_BH_tSNE: Principal component analysis')
			pca = decomposition.PCA(n_components=pca_dimensions)
			pca_matrix = pca.fit_transform(k_mer_frequency_matrix)
		else:
			logger.info('run_BH_tSNE: Principle component analysis step skipped')

		# 5. BH-tSNE
		logger.info('run_BH_tSNE: BH-tSNE')

		# Adjust perplexity according to the number of data points
		# Took logic from tsne source code
		if (len(k_mer_frequency_matrix) - 1) < (3 * perplexity)  :
			perplexity = (float(len(k_mer_frequency_matrix) - 1) / 3) - 1

		logger.info(str(len(k_mer_frequency_matrix)) + ' data points')
		logger.info(str(len(k_mer_frequency_matrix[0])) + ' dimensions')

		if (len(k_mer_frequency_matrix[0]) > pca_dimensions) and (do_pca == True):
			X = np.array(pca_matrix)
		else:
			X = np.array(k_mer_frequency_matrix)
		bh_tsne_matrix = bh_sne(X, d=2, perplexity=perplexity, theta=0.5)

		logger.info('run_BH_tSNE: Outputting file ' + output_filename)
		output = open(output_filename, 'w')
		# We will add bh_tsne_x and bh_tsne_y columns to the contig table
		contig_table = open(contig_table_path, 'r')
		contig_table_lines = contig_table.read().splitlines()

		for i, line in enumerate(contig_table_lines):
			if i == 0:
				new_header = line + '\tbh_tsne_x\tbh_tsne_y\n'
				output.write(new_header)
			else:
				line_list = line.split('\t')
				current_contig = line_list[0]
				if current_contig in sequence_names:
					contig_index = sequence_names[current_contig]
					bh_tsne_x = str(bh_tsne_matrix[contig_index][0])
					bh_tsne_y = str(bh_tsne_matrix[contig_index][1])
					output.write(line + '\t' + bh_tsne_x + '\t' + bh_tsne_y + '\n')

		output.close()

parser = argparse.ArgumentParser(description="Prototype script to automatically carry out secondary clustering of BH_tSNE coordinates based on DBSCAN and cluster purity")
parser.add_argument('-m','--marker_tab', help='Output of make_marker_table.py', required=True)
parser.add_argument('-d','--domain', help='Microbial domain (bacteria|archaea)', default='bacteria')
parser.add_argument('-f','--fasta', help='Assembly FASTA file', required=True) 
parser.add_argument('-o','--outdir', help='Path of directory for output', required=True)
parser.add_argument('-c','--contigtable', help='Output of make_contig_table.py', required=True)
parser.add_argument('-p','--processors', help='Number of processors used', default=1)
parser.add_argument('-t','--taxonomy_tab', help='Output of make_taxonomy_table.py')
args = vars(parser.parse_args())

hmm_table_path = args['marker_tab']
domain = args['domain']
fasta_path = args['fasta']
outdir = os.path.abspath(args['outdir'])
contig_table = args['contigtable']
processors = args['processors']
taxonomy_table_path = args['taxonomy_tab']

start_time = time.time()
FNULL = open(os.devnull, 'w')

pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)



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

# 2. If taxonomy data is available, load the information
taxonomy_info = dict()
if taxonomy_table_path:
	taxonomy_table = pd.read_table(taxonomy_table_path)
	for i,row in taxonomy_table.iterrows():
		taxonomy_info[row['contig']] = { 'phylum': row['phylum'], 'class': row['class'], 'order': row['order'], 'family': row['family'], 'genus': row['genus'], 'species': row['species']}

	# Make combined contig table
	new_contig_table_path = contig_table + '.taxonomy'
	new_contig_table = open(new_contig_table_path, 'w')
	with open(contig_table, 'r') as old_contig_table:
		for i,line in enumerate(old_contig_table):
			if i == 0:
				new_contig_table.write(line.rstrip() + '\tphylum\tclass\torder\tfamily\tgenus\tspecies\n')
			else:
				line_list = line.rstrip().split('\t')
				contig = line_list[0]
				new_contig_table.write(line.rstrip() + '\t' + taxonomy_info[contig]['phylum'] + '\t' + taxonomy_info[contig]['class'] + '\t' + taxonomy_info[contig]['order'] + '\t' + taxonomy_info[contig]['family'] + '\t' + taxonomy_info[contig]['genus'] + '\t' + taxonomy_info[contig]['species'] + '\n')
	new_contig_table.close()
	contig_table = new_contig_table_path

# Make a master table that will be updated as clustering is done

global_cluster_info = {}
global_cluster_contigs = {}
completeness_cutoff = 20
purity_cutoff = 100
round_counter = 0
current_fasta = fasta_path
BH_tSNE_counter = 0
master_table = None

# Load fasta sequences
assembly_seqs = {}
for seq_record in SeqIO.parse(fasta_path, 'fasta'):
	assembly_seqs[seq_record.id] = seq_record


# Run BH_tSNE
BH_tSNE_counter += 1
logger.info('Running BH-tSNE round ' + str(BH_tSNE_counter))
current_BH_tSNE_output = 'BH_tSNE' + str(BH_tSNE_counter) + '.tab'
# Carry out first BH_tSNE run
run_BH_tSNE(current_fasta,current_BH_tSNE_output,contig_table)

abs_BH_tSNE_path = os.path.abspath(current_BH_tSNE_output)
#BH_tSNE_r = get_table(abs_BH_tSNE_path)

if BH_tSNE_counter == 1:
	#master_table = pandas2ri.ri2py(BH_tSNE_r)
	master_table = pd.read_table(abs_BH_tSNE_path)

# Output current BH_tSNE table
#BH_tSNE_pd = pandas2ri.ri2py(BH_tSNE_r)
#BH_tSNE_pd.to_csv(path_or_buf=BH_tSNE_output_path, sep="\t", index=False, quoting=csv.QUOTE_NONE)

#current_r_table = BH_tSNE_r
current_table = pd.read_table(abs_BH_tSNE_path)



# Now if we have taxonomy data we do another round of clustering
if taxonomy_info:
	local_BH_tSNE_round = 0

	taxonomic_levels = ['phylum', 'class', 'order', 'family', 'genus', 'species']
	logger.info('Further splitting according to taxonomic classifications')
	for taxonomic_level in taxonomic_levels:
		logger.info('Taxonomic level: ' + taxonomic_level)
		combined_unclustered_table = pd.DataFrame()

		# Make subtables for each type of classification at the current level
		classification_dict = dict()
		for i,row in current_table.iterrows():
			classification_dict[row[taxonomic_level]] = 1

		# Skip iteration if the current taxonomic level is empty
		if not classification_dict:
			continue

		for classification in classification_dict.keys():
			logger.info('Examining ' + classification)
			# Get subset table
			subset_table = current_table.loc[current_table[taxonomic_level] == classification]

			while True:
				if len(subset_table.index) < 1:
					break
				round_counter += 1
				local_BH_tSNE_round += 1
				logger.info('Running DBSCAN round ' + str(round_counter))
				db_tables = runDBSCANs(subset_table)
				cluster_information, contig_cluster_dictionary, unclustered_table = assessDBSCAN(db_tables, contig_markers, domain, completeness_cutoff, purity_cutoff)

				subset_table = unclustered_table

				if not cluster_information:
					break

				# Populate the global data structures
				for	cluster in cluster_information:
					new_cluster_name = 'BH_tSNE' + str(BH_tSNE_counter) + '_round' + str(round_counter) + '_' + str(cluster)
					global_cluster_info[new_cluster_name] = cluster_information[cluster]

				for contig in contig_cluster_dictionary:
					new_cluster_name = 'BH_tSNE' + str(BH_tSNE_counter) + '_round' + str(round_counter) + '_' + str(contig_cluster_dictionary[contig])
					global_cluster_contigs[contig] = new_cluster_name

			# Add unclustered_table to combined unclustered dataframe
			combined_unclustered_table = combined_unclustered_table.append(unclustered_table)

		current_table = copy.deepcopy(combined_unclustered_table)
else:
	local_BH_tSNE_round = 0
	while True:
		round_counter += 1
		local_BH_tSNE_round += 1
		logger.info('Running DBSCAN round ' + str(round_counter))
	#db_tables = runDBSCANs(current_r_table)
		db_tables = runDBSCANs(current_table)
		cluster_information, contig_cluster_dictionary, unclustered_table = assessDBSCAN(db_tables, contig_markers, domain, completeness_cutoff, purity_cutoff)
		current_table = unclustered_table

		if not cluster_information:
			break

		# Populate the global data structures
		for	cluster in cluster_information:
			new_cluster_name = 'BH_tSNE' + str(BH_tSNE_counter) + '_round' + str(round_counter) + '_' + str(cluster)
			global_cluster_info[new_cluster_name] = cluster_information[cluster]

		for contig in contig_cluster_dictionary:
			new_cluster_name = 'BH_tSNE' + str(BH_tSNE_counter) + '_round' + str(round_counter) + '_' + str(contig_cluster_dictionary[contig])
			global_cluster_contigs[contig] = new_cluster_name

# If we are not done, write a new fasta for the next BH_tSNE run
# First convert the r table to a pandas table
#unclustered_pd = pandas2ri.ri2py(current_r_table)
unclustered_seqrecords = []
for i, row in unclustered_table.iterrows():
	contig = row['contig']
	unclustered_seqrecords.append(assembly_seqs[contig])

current_fasta = 'unclustered_for_BH_tSNE' + str(BH_tSNE_counter + 1) + '.fasta'
SeqIO.write(unclustered_seqrecords, current_fasta, 'fasta')


# Add cluster to master data frame
clusters = []
for contig in master_table['contig']:
	if contig in global_cluster_contigs:
		clusters.append(global_cluster_contigs[contig])
	else:
		clusters.append('unclustered')

master_table['cluster'] = clusters

# Output master_table
master_table_path = outdir + '/full_table'
master_table.to_csv(path_or_buf=master_table_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)

# Make summary table
summary_info = getClusterSummaryInfo(master_table)

# print summary table
summary_table_path = outdir + '/summary_table'
summary_table = open(summary_table_path, 'w')
summary_table.write('cluster\tsize\tlongest_contig\tn50\tnumber_contigs\tcompleteness\tpurity\tcov\tgc_percent\n')
for cluster in summary_info:
	size = str(summary_info[cluster]['size'])
	longest_contig = str(summary_info[cluster]['longest_contig'])
	n50 = str(summary_info[cluster]['n50'])
	number_contigs = str(summary_info[cluster]['number_contigs'])
	if cluster in global_cluster_info:
		completeness = str(global_cluster_info[cluster]['completeness'])
		purity = str(global_cluster_info[cluster]['purity'])
	else:
		completeness = 'NA'
		purity = 'NA'
	cov = str(summary_info[cluster]['cov'])
	gc = str(summary_info[cluster]['gc_percent'])
	status = None

	outputString = '\t'.join([str(cluster), size, longest_contig, n50, number_contigs, completeness, purity, cov, gc]) + '\n'

	summary_table.write(outputString)
summary_table.close()

# Now for each 'complete' cluster, we output a separate fasta file.  We collect all contigs from the 'failed' and 'noise' clusters and 
# output them to one fasta (unclustered.fasta)

# Initialize output fasta lists
output_fastas = {} # Keyed by cluster

for i, row in master_table.iterrows():
	contig = row['contig']
	cluster = row['cluster']
	seq_record = assembly_seqs[contig]
	if cluster in output_fastas:
		output_fastas[cluster].append(seq_record)
	else: 
		output_fastas[cluster] = [ seq_record ]

# Now write fasta files
for cluster in output_fastas:
	if cluster == 'unclustered':
		fasta_output_path = outdir + '/' + str(cluster) + '.fasta'
	else:
		fasta_output_path = outdir + '/cluster_' + str(cluster) + '.fasta'
	SeqIO.write(output_fastas[cluster], fasta_output_path, 'fasta')









