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
import random
#For multiprocessing

parser = argparse.ArgumentParser(description="Recruit unclustered (or non-marker)\
    sequences with Machine Learning classification using clustered sequence\
    as training data. Features to train with include sequence coverage,\
    composition, and homology. Confidence is calculated using jackknife\
    cross-validation by randomly subsetting the training data n number of times.")
parser.add_argument('-t','--contig_tab', help='Master contig table', required=True)
parser.add_argument('-c','--cluster_column', help='Name of column for cluster', \
    default='cluster')
#parser.add_argument('-r','--recursive', help='If specified, will run classification \
#    iteratively and refine traning data after each iteration.', action='store_true')
parser.add_argument('-C','--Confidence_cutoff', help='Confidence cutoff value\
    to use to keep ML-based predictions.', default=100)
parser.add_argument('-u','--unclustered_name', help='Name of unclustered group \
    in cluster column', default="unclustered")
parser.add_argument('-n','--num_iterations', help='Number of iterations for \
    jackknife cross-validation.', default=10)
parser.add_argument('-m','--k_mer_matrix', help='k-mer_matrix file.', default="k-mer_matrix")
parser.add_argument('-o','--out_table', help='Output table with new column\
    for ML-recruited sequences.',required=True)
args = vars(parser.parse_args())

args = {}
#In: /Users/Ian/Desktop/Bioinformatics_2017/autometa_testing/CHTC_parallelization/78.125Mbp/
args['contig_tab'] = "dbscan_2_3_dimensions.tab"
args['cluster_column'] = "cluster"
args['Confidence_cutoff'] = 100
args['k_mer_matrix'] = 'k-mer_matrix'
args['out_table'] = "test.out"
args['num_iterations'] = 10
args['unclustered_name'] = "unclustered"
args['unclustered_contig_list'] = "still_unclassified.list"

#############################
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
                print("-->This prediction adds marker redundancy...skipping...")
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
####################### 1.  Load contig table, kmer matrix, parse taxonomy info   ############################

#Mark start time
start_time = time.time()
#1. Load table with cluster info, taxonomy as binary matrices with pandas
print("Loading contig table...")
#Disable pandas warnings
pd.options.mode.chained_assignment = None
contig_table = pd.read_csv(args['contig_tab'],sep="\t")

print("Looking for taxonomy info in {}".format(args['contig_tab']))
use_taxonomy_info = False
try:
    phylum_dummy_matrix = pd.get_dummies(contig_table['phylum'])
    class_dummy_matrix = pd.get_dummies(contig_table['class'])
    order_dummy_matrix = pd.get_dummies(contig_table['order'])
    family_dummy_matrix = pd.get_dummies(contig_table['family'])
    genus_dummy_matrix = pd.get_dummies(contig_table['genus'])
    species_dummy_martix = pd.get_dummies(contig_table['species'])
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

bootstrap_iterations = int(args['num_iterations'])
confidence_cutoff = float(args['Confidence_cutoff'])
if confidence_cutoff % bootstrap_iterations != 0  and len(str(int(confidence_cutoff))) == len(str(bootstrap_iterations)):
    confidence_cutoff = round_down(confidence_cutoff,bootstrap_iterations)
cluster_column_name = args['cluster_column']
unclustered_name = args['unclustered_name']
taxonomy_matrix_dict = {}


####################### 2.  Load training data from contig table   ############################

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


####################### 3.  Train classifier and make predictions   ############################

iteration_start_time = time.time()
ML_predictions_dict = {}
ML_recruitment_list = []
recruited_sequence_length = 0
temp_contig_table = contig_table.copy(deep=True)
num_unclustered_contigs = contig_table[cluster_column_name].tolist().count(unclustered_name)
unclustered_contig_feature_list = []

#This is where the input list would be read in
if args['unclustered_contig_list'] != None:
    unclustered_contig_list_file = args['unclustered_contig_list']
    unclustered_contig_table = pd.DataFrame(columns = contig_table.columns)
    with open(unclustered_contig_list_file) as infile:
        for line in infile:
            contig = line.rstrip()
            global_contig_index = contig_index_dict[contig]
            unclustered_contig_table = unclustered_contig_table.append(contig_table.iloc[global_contig_index])
else:
    unclustered_contig_table = contig_table.loc[contig_table[cluster_column_name] == unclustered_name]
print("Recruiting {} unclustered sequences with {} training contigs. This could take a while...".format(len(unclustered_contig_table),len(features)))

newDF = pd.DataFrame(columns = contig_table.columns)

for count,row in unclustered_contig_table.iterrows():
    contig = row['contig']
    single_np_array = np.array([contig_feature_dict[contig]])
    contig_length = contig_table.iloc[count]['length']

    #Make predictions
    prediction,confidence = calculate_bootstrap_replicates(features, bootstrap_iterations)
    print("ML predictions and jackknife confidence for contig {}: {},{}".format(contig, prediction,confidence))
    redundant,is_marker_contig = redundant_marker_prediction(contig,prediction,temp_contig_table,cluster_column_name)

    if confidence >= confidence_cutoff and not redundant:
        #temp_contig_table[cluster_column_name].iloc[count] = prediction
        row[cluster_column_name] = prediction
        newDF = newDF.append(row)#, ignore_index = True)
    else:
        with open("still_unclassified.list", "a") as outfile:
            outfile.write(contig + "\n")

newDF.to_csv("new_classified.tab",sep="\t",index=False,mode="w",header=False,na_rep="NA")


"""
Output should be:
(1) pandas df for all classified contigs (if it doesn't exist already)
(2) pandas df for newly classified data (that could be concatenated to above df)
(3) list of still unclassified contigs
"""
