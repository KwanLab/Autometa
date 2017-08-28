#!/usr/bin/env python

import pandas as pd
import copy
import numpy as np
import csv
from sklearn import decomposition
from sklearn import tree,metrics,preprocessing
from sklearn.model_selection import train_test_split
import argparse
import sys
from joblib import Parallel, delayed
import pprint
import logging
import math
from scipy import stats
import collections
import pdb

def normalizeKmers(count_matrix): # list of lists, not a np matrix
	# We now remove all the k_mers where all counts are '1'
	print('Trimming k_mers')
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
	print('Normalizing counts')

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

def ML_assessClusters(table, confidence_cutoff = 50, singleton_cutoff = 90, clusters_to_examine = None):
	cluster_counts = dict() # Keyed by cluster, holds totals for 'congruent' classification or 'different' classification
	contig_reassignments = dict()

	subset_table = table.loc[table['cluster'] != 'unclustered']

	# Identify singletons
	singletons = findSingletons(table)

	# For each contig, we want to extract out its features, and make a custom training set that EXCLUDES that contig
	contig_list = list()
	cluster_list = list()
	ML_inputs = list()

	for i,row in subset_table.iterrows():
		current_contig = row['contig']
		current_cluster = row['cluster']

		if clusters_to_examine is not None:
			if current_cluster not in clusters_to_examine:
				continue

		contig_list.append(current_contig)
		cluster_list.append(current_cluster)

		if current_cluster not in cluster_counts:
			cluster_counts[current_cluster] = { 'congruent': 0, 'different': 0 }

	# Now do the ML prediction in parallel
	parallel_output = Parallel(n_jobs=processors)(delayed(reclassify_contig)(subset_table, contig_name, 10) for contig_name in contig_list)

	# Now we process the results
	for i,output_tuple in enumerate(parallel_output):
		ML_prediction,confidence = output_tuple
		current_cluster = cluster_list[i]
		current_contig = contig_list[i]

		if current_cluster in singletons:
			redundant = redundant_marker_prediction(current_contig, ML_prediction, subset_table, 'cluster')
			if confidence >= singleton_cutoff and not redundant:
				cluster_counts[current_cluster]['different'] += 1
				logger.debug('assessClusters: contig ' + current_contig + ', current cluster: ' + current_cluster + ', predicted: ' + ML_prediction + ', confidence: ' + str(confidence))
				contig_reassignments[current_contig] = 'unclustered'
			else:
				cluster_counts[current_cluster]['congruent'] += 1
		elif ML_prediction == current_cluster and confidence >= confidence_cutoff:
			cluster_counts[current_cluster]['congruent'] += 1
		else:
			cluster_counts[current_cluster]['different'] += 1
			logger.debug('assessClusters: contig ' + current_contig + ', current cluster: ' + current_cluster + ', predicted: ' + ML_prediction + ', confidence: ' + str(confidence))
			contig_reassignments[current_contig] = 'unclustered'

	# Determine congruent fractions
	cluster_results = dict()
	for cluster in cluster_counts:
		total_count = cluster_counts[cluster]['congruent'] + cluster_counts[cluster]['different']
		percentage = (float(cluster_counts[cluster]['congruent'])/total_count)*100
		cluster_results[cluster] = percentage

	return cluster_results, contig_reassignments

def reassignClusters(table, reassignments):
	# Now we go through the table and make the reassignments
	for i, row in table.iterrows():
		current_contig = row['contig']
		if current_contig in reassignments:
			if reassignments[current_contig] == 'unclustered':
				table.ix[i, 'cluster'] = reassignments[current_contig]
			else:
				redundant = redundant_marker_prediction(current_contig, reassignments[current_contig], table, 'cluster')
				logger.debug(current_contig + ' predicted to be in cluster ' + reassignments[current_contig] + ' but add redundancy')
				if not redundant:
					table.ix[i, 'cluster'] = reassignments[current_contig]

def findSingletons(table):
	cluster_contig_counts = dict()
	for i,row in table.iterrows():
		current_cluster = row['cluster']
		if current_cluster in cluster_contig_counts:
			cluster_contig_counts[current_cluster] += 1
		else:
			cluster_contig_counts[current_cluster] = 1

	singletons_found = dict()
	for cluster in cluster_contig_counts:
		if cluster_contig_counts[cluster] == 1:
			singletons_found[cluster] = 1

	return singletons_found

def reclassify_contig(table, contig_name, iterations = 10):
	# Set up data structures for training/classification
	to_be_classified_features = None

	features = list()
	labels = list()

	for j,row in table.iterrows():
		current_contig = row['contig']
		if has_taxonomy_info:
			current_features = np.array(taxonomy_matrix[j] + pca_matrix[j].tolist() + [ coverage_list[j] ])
		else:
			current_features = pca_matrix[j] + [ coverage_list[j] ]

		if current_contig == contig_name:
			to_be_classified_features = np.array([current_features])
		else:
			features.append(current_features)
			labels.append(row['cluster'])

	temp_dict = { 'training_features': features, 'training_labels': labels, 'classification_features': to_be_classified_features}
	classification, confidence = calculate_bootstrap_replicates(temp_dict, iterations)
	return classification, confidence

def redundant_marker_prediction(contig_name,predicted_cluster,pandas_table,cluster_column_name):
    #Function to check for redundancy of single copy gene markers in ML
    #predictions. The way it's written now,
    #Get a list of PFAM markers from current cluster
    cluster_df = pandas_table.loc[pandas_table[cluster_column_name] == predicted_cluster]
    cluster_PFAMs = []
    for count,contig in enumerate(cluster_df['contig']):
        #If it's a marker contig
        if cluster_df['num_single_copies'].iloc[count] > 0:
            contig_PFAMs = cluster_df['single_copy_PFAMs'].iloc[count].split(",")
            cluster_PFAMs += contig_PFAMs
    contig_index = pandas_table[pandas_table['contig'] == contig_name].index.tolist()
    #contig_df = pandas_table.loc[(pandas_table.contig == contig_name)]
    redundancy = False
    #if len(contig_PFAMs) > 0:
    #Avoid non-marker contigs in this calculation, for some reason evaluating to floats..
    if not isinstance(pandas_table['single_copy_PFAMs'][contig_index[0]],float):
        contig_PFAMs = pandas_table['single_copy_PFAMs'][contig_index[0]].split(",")
        for PFAM in contig_PFAMs:
            if PFAM in cluster_PFAMs:
                redundancy = True
                print("-->This prediction adds marker redundancy...skipping...")
                return redundancy
            else:
                pass
    else:
        #If no markers, can't add contamination...
        redundancy = False
        return redundancy
    #Maybe this update should happen after function return?
    #pandas_table[cluster_column_name][contig_index] = predicted_cluster
    return redundancy

def calculate_bootstrap_replicates(input_dictionary,iterations = 10):
	prediction_list = []
	training_features = input_dictionary['training_features']
	training_labels = input_dictionary['training_labels']
	classification_features = input_dictionary['classification_features']
	for i in range(iterations):
		#Here features and labels are global variables
		jackknifed_classifier = jackknife_training(training_features,training_labels)
		ML_prediction = jackknifed_classifier.predict(classification_features)[0]
		prediction_list.append(ML_prediction)
	counter = collections.Counter(prediction_list)
	top_prediction_set = counter.most_common(1)
	top_prediction = top_prediction_set[0][0]
	confidence = top_prediction_set[0][1]
	confidence_percent = round(float(confidence)/iterations*100,3)
	#To see frequency of all prediction: print counter
	return top_prediction,confidence_percent

def jackknife_training(features,labels):
	#Function to randomly subsample data into halves (hence 0.5), train
	#ML-classifier and make prediction. Used iteratively in
	#calculate_bootstap_replicates() function (see below)
	train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.50)
	my_classifier = tree.DecisionTreeClassifier()
	my_classifier = my_classifier.fit(train_features,train_labels)
	predictions = my_classifier.predict(test_features)
	return my_classifier

def revcomp( string ):
	trans_dict = { 'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C' }
	complement_list = list()

	for i in range(0, len(string) ):
		if string[i] in trans_dict:
			complement_list.append(trans_dict[string[i]])
		else:
			return -1

	return ''.join(reversed(complement_list))

parser = argparse.ArgumentParser(description="Perform initial clustering via BH-tSNE and DBSCAN.")
parser.add_argument('-t','--input_table', help='Master contig table. Optionally can contain taxonomy data', required=True)
parser.add_argument('-a','--assembly_fasta', help='Assembly fasta', required=True)
parser.add_argument('-o','--output_table', help='Path to output table', required=True)
parser.add_argument('-k','--kingdom', help='Kingdom to consider (archaea|bacteria)', choices=['bacteria','archaea'], default = 'bacteria')
parser.add_argument('-p','--processors', help='Number of processors used', default=1)
parser.add_argument('-m','--k_mer_matrix', help='Path to k_mer matrix file', required=True)

args = vars(parser.parse_args())

input_table_path = args['input_table']
fasta_path = args['assembly_fasta']
output_table_path = args['output_table']
domain = args['kingdom']
processors = int(args['processors'])
matrix_file = args['k_mer_matrix']

#logger
logger = logging.getLogger('ML_prune.py')
hdlr = logging.FileHandler('ML_prune.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
console = logging.StreamHandler()
console.setFormatter(formatter)
logger.setLevel(logging.DEBUG)
logger.addHandler(console)

pp = pprint.PrettyPrinter(indent=4)

pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)

master_table = pd.read_table(input_table_path)

# Copy single-copy marker information to a data structure
contig_markers = {}
for i, row in master_table.iterrows():
	pfamString = row['single_copy_PFAMs']
	contig = row['contig']
	if pfamString == 'NA' or (isinstance(pfamString, float) and math.isnan(pfamString)):
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

# Now we load the k_mer matrix
k_mer_dict = dict() # Holds lists of k_mer counts, keyed by contig name

k_mer_size = 5
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

with open(matrix_file) as matrix:
	for i,line in enumerate(matrix):
		if i > 0:
			line_list = line.rstrip().split('\t')
			contig = line_list.pop(0)
			line_list = [ int(x) for x in line_list ]
			k_mer_dict[contig] = line_list

### Collate training data for ML steps later
# We now set up global data structures to be used in supervised machine learning
contig_list = master_table['contig'].tolist()
coverage_list = master_table['cov'].tolist()

if 'kingdom' in master_table.columns:
	has_taxonomy_info = True
else:
	has_taxonomy_info = False

taxonomy_matrix = list()

if has_taxonomy_info:
	phylum_dummy_matrix = pd.get_dummies(master_table['phylum'])
	class_dummy_matrix = pd.get_dummies(master_table['class'])
	order_dummy_matrix = pd.get_dummies(master_table['order'])
	family_dummy_matrix = pd.get_dummies(master_table['family'])
	genus_dummy_matrix = pd.get_dummies(master_table['genus'])
	species_dummy_martix = pd.get_dummies(master_table['species'])

	for j,contig in enumerate(master_table['contig']):
		tax_phylum = list(phylum_dummy_matrix.iloc[j])
		tax_class = list(class_dummy_matrix.iloc[j])
		tax_order = list(order_dummy_matrix.iloc[j])
		tax_family = list(family_dummy_matrix.iloc[j])
		tax_genus = list(genus_dummy_matrix.iloc[j])
		tax_species = list(species_dummy_martix.iloc[j])
		taxonomy = tax_phylum + tax_class + tax_order + tax_family + tax_genus + tax_species
		taxonomy_matrix.append(taxonomy)

# Make normalized k_mer matrix
k_mer_counts = list()
for contig in contig_list:
	k_mer_counts.append(k_mer_dict[contig])

normalized_k_mer_matrix = normalizeKmers(k_mer_counts)

# For performance reasons we reduce the dimensions to 50 with PCA
pca = decomposition.PCA(n_components=50)
pca_matrix = pca.fit_transform(normalized_k_mer_matrix)

### Pruning of clusters through supervised machine learning classification
### Based on whether each contig re-classifies to the same cluster when it is taken out of the training set

# Now we refine the clustering using supervised ML
all_good_clusters = False
iteration = 0

# Do initial iteration, where we generate the first scores
print('Cluster pruning iteration: ' + str(iteration))
logger.info('Cluster pruning iteration: ' + str(iteration))

cluster_scores, contig_reassignments = ML_assessClusters(master_table) # cluster_scores pre-sorted in ascending order

logger.debug(pprint.pformat(cluster_scores))
iteration += 1

good_clusters = dict()
bad_clusters = True

while bad_clusters:
	clusters_to_consider = dict()
	# Make note of "good" clusters that are >= 100%
	for cluster in cluster_scores:
		if cluster_scores[cluster] == 100:
			good_clusters[cluster] = 1
		else:
			clusters_to_consider[cluster] = 1

	# Now filter out reassignments according to the running list of "good" clusters
	contig_clusters = dict()
	for i,row in master_table.iterrows():
		contig = row['contig']
		cluster = row['cluster']
		contig_clusters[contig] = cluster

	filtered_contig_reassignments = dict()
	for contig in contig_reassignments:
		current_cluster = contig_clusters[contig]
		if current_cluster in good_clusters:
			continue
		else:
			filtered_contig_reassignments[contig] = contig_reassignments[contig]

	reassignClusters(master_table, filtered_contig_reassignments)

	print('Cluster pruning iteration: ' + str(iteration))
	logger.info('Cluster pruning iteration: ' + str(iteration))

	cluster_scores, contig_reassignments = ML_assessClusters(master_table, 50, 90, clusters_to_consider)
	logger.debug(pprint.pformat(cluster_scores))
	iteration += 1

	bad_clusters = False
	for cluster in cluster_scores:
		if cluster_scores[cluster] < 90 and cluster not in good_clusters:
			bad_clusters = True


# Output table
master_table.to_csv(path_or_buf=output_table_path, sep='\t', index=False, quoting=csv.QUOTE_NONE)
