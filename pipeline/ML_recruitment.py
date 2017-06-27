#!/usr/bin/env python

from __future__ import division
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

parser = argparse.ArgumentParser(description="Recruit unclustered (or non-marker)\
    sequences with Machine Learning classification using clustered sequence\
    as training data. Features to train with include sequence coverage,\
    composition, and homology. Confidence is calculated using jackknife\
    cross-validation by randomly subsetting the training data n number of times.")
parser.add_argument('-t','--contig_tab', help='Master contig table', required=True)
parser.add_argument('-c','--cluster_column', help='Name of column for cluster', \
    default='cluster')
parser.add_argument('-C','--Confidence_cutoff', help='Confidence cutoff value\
    to use to keep ML-based predictions.', default=95)
parser.add_argument('-u','--unclustered_name', help='Name of unclustered group \
    in cluster column', default="unclustered")
parser.add_argument('-n','--num_iterations', help='Number of iterations for \
    jackknife cross-validation.', default=10)
parser.add_argument('-m','--k_mer_matrix', help='k-mer_matrix file.', default="k-mer_matrix")
parser.add_argument('-o','--out_table', help='Output table with new column\
    for ML-recruited sequences.',required=True)
args = vars(parser.parse_args())

#1. Load table with cluster info, taxonomy as binary matrices with pandas
print("Loading contig table..")
#Disable pandas warnings
pd.options.mode.chained_assignment = None
contig_table = pd.read_csv(args['contig_tab'],sep="\t")

print("Loading taxonomy info as dummy matrices..")
phylum_dummy_matrix = pd.get_dummies(contig_table['phylum'])
class_dummy_matrix = pd.get_dummies(contig_table['class'])
order_dummy_matrix = pd.get_dummies(contig_table['order'])
family_dummy_matrix = pd.get_dummies(contig_table['family'])
genus_dummy_matrix = pd.get_dummies(contig_table['genus'])
species_dummy_martix = pd.get_dummies(contig_table['species'])

def round_down(num, divisor):
    return num - (num%divisor)

###For k-kmer matrix reduction - START

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
#length_weight = []
for count,contig in enumerate(contig_table['contig']):
    #Actually the label here should be the bin name
    #known_genome =  contig_table['known_genome'][count]
    cluster = contig_table[cluster_column_name][count]
    tax_phylum = list(phylum_dummy_matrix.iloc[count])
    tax_class = list(class_dummy_matrix.iloc[count])
    tax_order = list(order_dummy_matrix.iloc[count])
    tax_family = list(family_dummy_matrix.iloc[count])
    tax_genus = list(genus_dummy_matrix.iloc[count])
    tax_species = list(species_dummy_martix.iloc[count])
    taxonomy = tax_phylum + tax_class + tax_order + tax_family + tax_genus + tax_species
    taxonomy_matrix_dict[contig] = taxonomy
    if cluster != unclustered_name:
        bh_tsne_x = contig_table['bh_tsne_x'][count]
        bh_tsne_y = contig_table['bh_tsne_y'][count]
        length = contig_table['length'][count]
        cov = contig_table['cov'][count]
        gc = contig_table['gc'][count]
        #Actually the label here should be the bin name
        #label = known_genome
        label = cluster
        #features.append([bh_tsne_x,bh_tsne_y,cov] + taxonomy)
        features.append(pca_matrix[count].tolist() + [cov] + taxonomy)
        labels.append(label)

def jackknife_training(features,labels):
    #Function to randomly subsample data into halves (hence 0.5), train
    #ML-classifier and make prediction. Used iteratively in
    #calculate_bootstap_replicates() function (see below)
    train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.50)
    my_classifier = tree.DecisionTreeClassifier()
    my_classifier = my_classifier.fit(train_features,train_labels)
    predictions = my_classifier.predict(test_features)
    return my_classifier

def calculate_bootstap_replicates(feature_array,iterations = 10):
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
    #if len(contig_PFAMs) > 0:
    #Avoid non-marker contigs in this calculation, for some reason evaluating to floats..
    if not isinstance(pandas_table.iloc[contig_index]['single_copy_PFAMs'],float):
        contig_PFAMs = pandas_table['single_copy_PFAMs'][contig_index].split(",")
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

num_confident_predictions = 1
iteration = 0
while num_confident_predictions > 0:
    ML_predictions_dict = {}
    ML_recruitment_list = []
    prediction_accuracy_list = []
    temp_contig_table = contig_table.copy(deep=True)
    #Recruit unclustered sequences
    num_unclustered_contigs = contig_table[cluster_column_name].tolist().count(unclustered_name)
    print("Recruiting {} unclustered sequences. This could take a while...".format(num_unclustered_contigs))
    for count,contig in enumerate(contig_table['contig']):
        taxonomy = taxonomy_matrix_dict[contig]
        bh_tsne_x = contig_table.iloc[count]['bh_tsne_x']
        bh_tsne_y = contig_table.iloc[count]['bh_tsne_y']
        composition_feature_list = pca_matrix[count].tolist() + [cov]
        #Concatenate composition and taxonomy features into single list
        single_np_array = np.array([composition_feature_list + taxonomy])
        contig_length = contig_table.iloc[count]['length']
        #After the first iteration, train from previous confident predictions
        if iteration >= 1:
            cluster = contig_table.iloc[count]['ML_expanded_clustering']
        #Otherwise, just use the initial clustering as labels
        else:
            cluster = contig_table.iloc[count][cluster_column_name]

        #If contig is unclustered, make ML prediction
        if cluster == unclustered_name:# and contig_table.iloc[count]['num_single_copies'] > 0:
            ML_prediction,confidence = calculate_bootstap_replicates(single_np_array,bootstrap_iterations)
            ML_predictions_dict[contig] = ML_prediction,confidence
            print("ML predictions and jackknife confidence for contig {}: {},{}".format(contig, ML_prediction,confidence))
            #If it the prediction passes confidence cutoff
            #Could also look for redundant markers...
            redundant = redundant_marker_prediction(contig,ML_prediction,temp_contig_table,cluster_column_name)
            if confidence >= confidence_cutoff and not redundant:
                #Add prediction to ML_recruitment_list
                ML_recruitment_list.append(ML_prediction)
                prediction_accuracy_list.append(ML_prediction)
                #Update contig table, so that any markers added to the cluster will
                #be considered in the next check of marker redundancy
                temp_contig_table[cluster_column_name].iloc[count] = ML_prediction
            else:
                ML_recruitment_list.append(unclustered_name)
        else:
            ML_recruitment_list.append(cluster)

    num_predictions = len(ML_predictions_dict)
    num_confident_predictions = len(prediction_accuracy_list)
    print("{} of {} predictions were {}% confident and non-redundant for iteration {}".format(num_confident_predictions,num_predictions,confidence_cutoff,iteration))

    contig_table['ML_expanded_clustering'] = ML_recruitment_list
    iteration += 1

#Write out final table
contig_table.to_csv(args['out_table'],sep="\t",index=False)
