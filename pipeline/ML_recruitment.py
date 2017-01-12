#!/usr/bin/env python

from __future__ import division
import numpy as np
import pandas as pd
from sklearn import tree,metrics,preprocessing
from sklearn.model_selection import train_test_split
import collections
import argparse

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
parser.add_argument('-m','--train_w_markers', help='Only train the ML classifier\
    with marker contigs.', default="True")
parser.add_argument('-o','--out_table', help='Output table with new column\
    for ML-recruited sequences.',required=True)
args = vars(parser.parse_args())

#Set load paramters - convert to argparse
bootstrap_iterations = int(args['num_iterations'])
confidence_cutoff = float(args['Confidence_cutoff'])
cluster_column_name = args['cluster_column']
unclustered_name = args['unclustered_name']
train_w_markers = args['train_w_markers']
taxonomy_matrix_dict = {}

#1. Load table with cluster info, taxonomy as binary matrices with pandas
print("Loading contig table..")
contig_table = pd.read_csv(args['contig_tab'],sep="\t")

print("Loading taxonomy info as dummy matrices..")
phylum_dummy_matrix = pd.get_dummies(contig_table['phylum'])
class_dummy_matrix = pd.get_dummies(contig_table['class'])
order_dummy_matrix = pd.get_dummies(contig_table['order'])
family_dummy_matrix = pd.get_dummies(contig_table['family'])
genus_dummy_matrix = pd.get_dummies(contig_table['genus'])
species_dummy_martix = pd.get_dummies(contig_table['species'])

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
    if train_w_markers.lower() == "true":
        if contig_table['num_single_copies'][count] > 0 and cluster != unclustered_name:
            vizbin_x = contig_table['vizbin_x'][count]
            vizbin_y = contig_table['vizbin_y'][count]
            length = contig_table['length'][count]
            cov = contig_table['cov'][count]
            gc = contig_table['gc'][count]
            #Actually the label here should be the bin name
            #label = known_genome
            label = cluster
            features.append([vizbin_x,vizbin_y,cov] + taxonomy)
            labels.append(label)
    elif train_w_markers.lower() == "false":
        if cluster != unclustered_name:
            vizbin_x = contig_table['vizbin_x'][count]
            vizbin_y = contig_table['vizbin_y'][count]
            length = contig_table['length'][count]
            cov = contig_table['cov'][count]
            gc = contig_table['gc'][count]
            #Actually the label here should be the bin name
            #label = known_genome
            label = cluster
            features.append([vizbin_x,vizbin_y,cov] + taxonomy)
            labels.append(label)
    else:
        print("Unrecognized argument for -m.\nExiting...")
        exit()

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

ML_predictions_dict = {}
ML_recruitment_list = []
prediction_accuracy_list = []
temp_contig_table = contig_table.copy(deep=True)
#Recruit unclustered sequences
print("Recruiting unclustered sequences. This could take a while...")
for count,contig in enumerate(contig_table['contig']):
    taxonomy = taxonomy_matrix_dict[contig]
    vizbin_x_x = contig_table.iloc[count]['vizbin_x']
    vizbin_y_x = contig_table.iloc[count]['vizbin_y']
    cov_x = contig_table.iloc[count]['cov']
    composition_feature_list = [vizbin_x_x,vizbin_y_x,cov_x]
    #Concatenate composition and taxonomy features into single list
    single_np_array = np.array([composition_feature_list + taxonomy])
    contig_length = contig_table.iloc[count]['length']
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
print("{} of {} predictions were {}% confident and non-redundant".format(num_confident_predictions,num_predictions,confidence_cutoff))

contig_table['ML_expanded_clustering'] = ML_recruitment_list
contig_table.to_csv(args['out_table'],sep="\t",index=False)
