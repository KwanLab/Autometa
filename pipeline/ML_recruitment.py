#!/usr/bin/env python

from __future__ import division
import time
import numpy as np
import pandas as pd
from sklearn import tree,metrics,preprocessing
from sklearn.model_selection import train_test_split
import collections
import argparse
#For kmer matrix reduction
from scipy import stats
import math
from sklearn import decomposition
#for parallel ML
from joblib import Parallel, delayed
import random
import multiprocessing
import os

parser = argparse.ArgumentParser(description="Recruit unclustered (or non-marker)\
    sequences with Machine Learning classification using clustered sequence\
    as training data. Features to train with include sequence coverage,\
    composition, and homology. Confidence is calculated using jackknife\
    cross-validation by randomly subsetting the training data n number of times.")
parser.add_argument('-t','--contig_tab', metavar='<contig.tab>', help='Path to master contig table which includes initial clusters', required=True)
parser.add_argument('-c','--cluster_column', metavar='<column header>', help='Name of column containing initial cluster information', \
    default='cluster')
parser.add_argument('-p','--processors', metavar='<int>', help='Number of processors to use', type=int, \
    default=1)
parser.add_argument('-r','--recursive', help='If specified, will run classification \
    iteratively and refine traning data after each iteration.', action='store_true')
parser.add_argument('-C','--Confidence_cutoff', metavar='<int>', help='Confidence cutoff value\
    to use to keep ML-based predictions.', type=int, default=100)
parser.add_argument('-u','--unclustered_name', metavar='<unclustered name>', help='Name of unclustered group \
    in cluster column', default="unclustered")
parser.add_argument('-n','--num_iterations', metavar='<int>', help='Number of iterations for \
    jackknife cross-validation.', type=int, default=10)
parser.add_argument('-m','--k_mer_matrix', metavar='<k-mer.tab>', help='Path to k-mer_matrix file.', default="k-mer_matrix")
parser.add_argument('-o','--out_table', metavar='<output.tab>', help='Path to create output table with new column\
    for ML-recruited sequences.',required=True)
parser.add_argument('-k','--kingdom', metavar='<archaea|bacteria>', help='Kingdom to consider (archaea|bacteria)',\
    choices=['bacteria','archaea'], default = 'bacteria')
args = vars(parser.parse_args())

def round_down(num, divisor):
    return num - (num%divisor)

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
	#logger.info('Trimming k-mers')
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
		#unique_k_mers is a global variable
		for j in range(0, len(unique_k_mers.keys())):
			if j not in columns_to_delete:
				new_row.append(count_matrix[i][j])
		filtered_contig_k_mer_counts.append(new_row)

	# 3. Normalization
	#logger.info('Normalizing counts')

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

def jackknife_training(features,labels):
    #Function to randomly subsample data into halves (hence 0.5), train
    #ML-classifier and make prediction. Used iteratively in
    #calculate_bootstrap_replicates() function (see below)
    train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.50)
    my_classifier = tree.DecisionTreeClassifier()
    my_classifier = my_classifier.fit(train_features,train_labels)
    predictions = my_classifier.predict(test_features)
    return my_classifier

def calculate_bootstrap_replicates(feature_array,iterations = 10):
    prediction_list = []
    for i in range(iterations):
        #Here features and labels are global variables
        jackknifed_classifier = jackknife_training(features,labels)
        ML_prediction = jackknifed_classifier.predict(feature_array)[0]
        prediction_list.append(ML_prediction)
    counter = collections.Counter(prediction_list)
    top_prediction_set = counter.most_common(1)
    top_prediction = top_prediction_set[0][0]
    confidence = top_prediction_set[0][1]
    confidence_percent = round(confidence/iterations*100,3)
    #To see frequency of all prediction: print counter
    return top_prediction,confidence_percent

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
    contig_index = list(pandas_table['contig']).index(contig_name)
    #contig_df = pandas_table.loc[(pandas_table.contig == contig_name)]
    redundancy = False
    is_marker_contig = False
    #Avoid non-marker contigs in this calculation, for some reason evaluating to floats..
    if not isinstance(pandas_table.iloc[contig_index]['single_copy_PFAMs'],float):
        contig_PFAMs = pandas_table['single_copy_PFAMs'][contig_index].split(",")
        if len(contig_PFAMs) > 0:
            is_marker_contig = True
        for PFAM in contig_PFAMs:
            if PFAM in cluster_PFAMs:
                redundancy = True
                #print("-->This prediction adds marker redundancy...skipping...")
                return redundancy,is_marker_contig
            else:
                pass
    else:
        #If no markers, can't add contamination...
        redundancy = False
        return redundancy,is_marker_contig
    #Maybe this update should happen after function return?
    #pandas_table[cluster_column_name][contig_index] = predicted_cluster
    return redundancy,is_marker_contig

def calculateClusterStats(pandas_table,cluster_column,life_domain="bacteria"):
    #Function to calculate completeness and contamination of every dbscan
    #clusters for a given clustering condition (e.g. eps value)
    if life_domain=="bacteria":
        expected_number = 139
    elif life_domain=="archaea":
        expected_number = 164
    else:
        print("Unexpected life domain: {}. Please select 'bacteria' or 'archaea'.".format(life_domain))
        exit()

    cluster_dict = {}
    #Will probably have to change some of this to accomodate non-marker contigs
    for count,contig in enumerate(pandas_table['contig']):
        label = pandas_table[cluster_column][count]
        num_markers = pandas_table['num_single_copies'][count]
        if num_markers >= 1:
            if label not in cluster_dict:
                cluster_dict[label] = {}
                cluster_dict[label]['PFAMs'] = pandas_table.iloc[count]['single_copy_PFAMs'].split(",")
            else:
                cluster_dict[label]['PFAMs'] += pandas_table.iloc[count]['single_copy_PFAMs'].split(",")

    #Now calculate completeness and purity for each cluster
    for cluster,PFAM_dict in cluster_dict.items():
        PFAM_list = PFAM_dict.values()[0]
        counter = collections.Counter(PFAM_list)
        completeness = len(counter)/expected_number*100

        repeated_markers = 0
        for PFAM,frequency in dict(counter).items():
            if frequency > 1:
                repeated_markers += 1

        purity = 100 - (repeated_markers/expected_number*100)
        cluster_dict[cluster]['completeness'] = completeness
        cluster_dict[cluster]['purity'] = purity

    return cluster_dict

#Mark start time
start_time = time.time()

# Do check to see if the contig table and the k_mer_matrix file exists
if not os.path.isfile(args['contig_tab']):
    print('Error! Could not find contig table at the following path: ' + args['contig_tab'])
    exit(1)

if not os.path.isfile(args['k_mer_matrix']):
    print('Error! Could not find k-mer matrix file at the following path: ' + args['k_mer_matrix'])
    exit(1)

#1. Load table with cluster info, taxonomy as binary matrices with pandas
print("Loading contig table...")
#Disable pandas warnings
pd.options.mode.chained_assignment = None
contig_table = pd.read_csv(args['contig_tab'],sep="\t")
kingdom = args['kingdom']

print("Looking for taxonomy info in {}".format(args['contig_tab']))
use_taxonomy_info = False
try:
    phylum_dummy_matrix = pd.get_dummies(contig_table['phylum']).astype(np.int8)
    class_dummy_matrix = pd.get_dummies(contig_table['class']).astype(np.int8)
    order_dummy_matrix = pd.get_dummies(contig_table['order']).astype(np.int8)
    family_dummy_matrix = pd.get_dummies(contig_table['family']).astype(np.int8)
    genus_dummy_matrix = pd.get_dummies(contig_table['genus']).astype(np.int8)
    species_dummy_martix = pd.get_dummies(contig_table['species']).astype(np.int8)
    print("Loaded taxonomy info as dummy matrices...")
    use_taxonomy_info = True
except KeyError:
    print("Couldn't find taxonomy info in table. Excluding as training feature...")

#contig_table = "full_table"
master_table = contig_table

#Define "unique_k_mers"
# Count K-mer frequencies
k_mer_size = 5
matrix_file = args['k_mer_matrix']
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

# Now we load the k-mer matrix
print("Loading k-mer matrix...")
k_mer_dict = {}
with open(matrix_file) as matrix:
	for i,line in enumerate(matrix):
		if i > 0:
			line_list = line.rstrip().split('\t')
			contig = line_list.pop(0)
			line_list = [ int(x) for x in line_list ]
			k_mer_dict[contig] = line_list

# Make normalized k-mer matrix
print("Normalizing k-mer martix...")
contig_list = master_table['contig'].tolist()
k_mer_counts = list()
for contig in contig_list:
	k_mer_counts.append(k_mer_dict[contig])

normalized_k_mer_matrix = normalizeKmers(k_mer_counts)

print("Reducing normalized k-mer matrix to 50 dimensions with PCA...")
# For performance reasons we reduce the dimensions to 50 with PCA
pca = decomposition.PCA(n_components=50)
pca_matrix = pca.fit_transform(normalized_k_mer_matrix)

###For k-kmer matrix reduction - END

#Set load paramters - convert to argparse
processors = int(args['processors'])
if processors > multiprocessing.cpu_count():
    processors = multiprocessing.cpu_count()
bootstrap_iterations = int(args['num_iterations'])
confidence_cutoff = float(args['Confidence_cutoff'])
if confidence_cutoff % bootstrap_iterations != 0  and len(str(int(confidence_cutoff))) == len(str(bootstrap_iterations)):
    confidence_cutoff = round_down(confidence_cutoff,bootstrap_iterations)
cluster_column_name = args['cluster_column']
unclustered_name = args['unclustered_name']
taxonomy_matrix_dict = {}

#2. Parse vizbin, cov, and taxonomy info in "features" and autometa-defined
# clusters into "labels" for classifier using appropriate data structure
print("Loading other features and labels...")
features = []
labels = []
contig_index_dict = {}
contig_feature_dict = {}
for count,contig in enumerate(contig_table['contig']):
    num_markers = contig_table['num_single_copies'][count]
    num_single_copies = contig_table
    contig_index_dict[contig] = count
    bh_tsne_x = contig_table['bh_tsne_x'][count]
    bh_tsne_y = contig_table['bh_tsne_y'][count]
    length = contig_table['length'][count]
    cov = contig_table['cov'][count]
    gc = contig_table['gc'][count]
    cluster = contig_table[cluster_column_name][count]
    contig_feature_dict[contig] = pca_matrix[count].tolist() + [cov]
    if use_taxonomy_info:
        tax_phylum = list(phylum_dummy_matrix.iloc[count])
        tax_class = list(class_dummy_matrix.iloc[count])
        tax_order = list(order_dummy_matrix.iloc[count])
        tax_family = list(family_dummy_matrix.iloc[count])
        tax_genus = list(genus_dummy_matrix.iloc[count])
        tax_species = list(species_dummy_martix.iloc[count])
        taxonomy = tax_phylum + tax_class + tax_order + tax_family + tax_genus + tax_species
        taxonomy_matrix_dict[contig] = taxonomy
        contig_feature_dict[contig] = pca_matrix[count].tolist() + [cov] + taxonomy
    if cluster != unclustered_name and num_markers > 0:
        features.append(contig_feature_dict[contig])
        labels.append(cluster)

print("There are {} training contigs...".format(len(features)))

num_confident_predictions = 1
num_markers_classifed = 1
iteration = 0
while num_markers_classifed > 0:
    classified_marker_list = []
    iteration_start_time = time.time()
    ML_predictions_dict = {}
    ML_recruitment_list = []
    recruited_sequence_length = 0
    accurate_prediction_list = []
    temp_contig_table = contig_table.copy(deep=True)
    #Recruit unclustered sequences
    if iteration > 0:
        cluster_column_name = "ML_expanded_clustering"
    num_unclustered_contigs = contig_table[cluster_column_name].tolist().count(unclustered_name)
    unclustered_contig_feature_list = []
    unclustered_contig_list = []
    print("Recruiting {} unclustered sequences with {} training contigs. This could take a while...".format(num_unclustered_contigs,len(features)))

    #Prepare unclustered contig feature array
    for count,contig in enumerate(contig_table['contig']):
        single_np_array = np.array([contig_feature_dict[contig]])
        contig_length = contig_table.iloc[count]['length']
        #After the first iteration, train from previous confident predictions
        cluster = contig_table.iloc[count][cluster_column_name]

        #If contig is unclustered, prepare feature array for (multiproccesed) ML prediction
        if cluster == unclustered_name:
            unclustered_contig_feature_list.append(single_np_array)
            unclustered_contig_list.append(contig)
            ML_recruitment_list.append(unclustered_name)

        #I think the index of the list is getting messed up by this, actually. Maybe store index and cluster assingment in separate df
        else:
            ML_recruitment_list.append(cluster)

    #START - Multiprocess subroutine
    multiprocessed_output = Parallel(n_jobs = processors)(delayed(calculate_bootstrap_replicates)(features, bootstrap_iterations) for features in unclustered_contig_feature_list)
    #END - Multiprocess subroutine
    for count,output_tuple in enumerate(multiprocessed_output):
        #Assuming index (from multiprocessing) is preserved
        contig = unclustered_contig_list[count]
        ML_prediction,confidence = output_tuple
        ML_predictions_dict[contig] = ML_prediction,confidence
        #print("ML predictions and jackknife confidence for contig {}: {},{}".format(contig, ML_prediction,confidence))
        #If it the prediction passes confidence cutoff
        #Could also look for redundant markers...
        redundant,is_marker_contig = redundant_marker_prediction(contig,ML_prediction,temp_contig_table,cluster_column_name)
        global_contig_index = contig_index_dict[contig]
        if confidence >= confidence_cutoff and not redundant:
            print("ML predictions and jackknife confidence for contig {}: {},{}".format(contig, ML_prediction,confidence))
            #Add prediction to ML_recruitment_list/replace with updated label
            #ML_recruitment_list.append(ML_prediction)
            ML_recruitment_list[global_contig_index] = ML_prediction
            accurate_prediction_list.append(ML_prediction)
            recruited_sequence_length += contig_length
            #Update contig table, so that any markers added to the cluster will
            #be considered in the next check of marker redundancy
            temp_contig_table[cluster_column_name].iloc[count] = ML_prediction
            contig_table[cluster_column_name].iloc[global_contig_index] = ML_prediction #NOTE: Think this may be source of error

            #Update training data with any confident and non-redundant marker contig classifications
            if is_marker_contig:
                features.append(contig_feature_dict[contig])
                labels.append(ML_prediction)
                classified_marker_list.append(ML_prediction)
        else:
            #ML_recruitment_list.append(unclustered_name)
            ML_recruitment_list[global_contig_index] = unclustered_name

    num_predictions = len(multiprocessed_output)
    num_confident_predictions = len(accurate_prediction_list)
    num_markers_classifed = len(classified_marker_list)
    #Calculate average cluster stats
    cluster_stats_dict = calculateClusterStats(contig_table,cluster_column_name,kingdom)
    completeness_list = []
    purity_list = []
    for cluster,info_dictionary in cluster_stats_dict.items():
        completeness_list.append(info_dictionary['completeness'])
        purity_list.append(info_dictionary['purity'])
    mean_completeness = round(np.mean(completeness_list),1)
    mean_purity = round(np.mean(purity_list),1)
    elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - iteration_start_time),2)))
    print("{} ({} marker contigs) of {} predictions ({} bp) were {}% confident and non-redundant for iteration {} in {} (HH:MM:SS). Mean completeness,purity: {},{}"\
        .format(num_confident_predictions,num_markers_classifed,num_predictions,recruited_sequence_length,confidence_cutoff,iteration,elapsed_time,mean_completeness,mean_purity))

    contig_table['ML_expanded_clustering'] = ML_recruitment_list
    #Break after first iteration if recursive option not specified:
    if not args['recursive']:
        break
    iteration += 1

elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))
print("Done! Total elapsed time = {} (HH:MM:SS)".format(elapsed_time))
#Write out final table
contig_table.to_csv(args['out_table'],sep="\t",index=False)
