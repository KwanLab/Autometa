#!/usr/bin/env python

from __future__ import division
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn import tree,metrics,preprocessing
from sklearn.cross_validation import train_test_split
from sklearn.metrics import accuracy_score
import collections
import argparse

parser = argparse.ArgumentParser(description="Recruit unclustered or non-marker\
    sequences with Machine Learning classification using clustered sequence\
    as training data. Features to train with include sequence coverage,\
    composition, and homology. Confidence is calculated using jackknife\
    cross-validation by randomly subsetting the data X number of times.")
parser.add_argument('-t','--contig_tab', help='Contig table', required=True)
parser.add_argument('-c','--cluster_column', help='Name of column for cluster', \
    default='cluster')
parser.add_argument('-C','--Confidence_cutoff', help='Confidence cutoff value\
    to use to keep ML-based predictions.', default=95)
parser.add_argument('-u','--unclustered_name', help='Name of unclustered group \
    in cluster column', default="unclustered")
parser.add_argument('-X','--num_iterations', help='Number of iterations for \
    jackknife cross-validation.', default=100)
parser.add_argument('-o','--out_table', help='Output table with new column\
    for ML-recruited sequences.',required=True)
args = vars(parser.parse_args())

#Set load paramters - convert to argparse
bootstrap_iterations = args['num_iterations']
confidence_cutoff = args['Confidence_cutoff']
cluster_column_name = args['cluster_column']
unclustered_name = args['unclustered_name']
taxonomy_matrix_dict = {}

#1. Load table with cluster info, taxonomy as binary matrices with pandas
print("Loading contig table..")
contig_table = pd.read_csv("../MIX_51_master_w_taxonomy.tab",sep="\t")

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
short_features = [] #Short classifier features
labels = []
length_weight = []
for count,contig in enumerate(contig_table['contig']):
    #Only train with marker contigs from known genome bins
    known_genome =  contig_table['known_genome'][count]
    tax_phylum = list(phylum_dummy_matrix.iloc[count])
    tax_class = list(class_dummy_matrix.iloc[count])
    tax_order = list(order_dummy_matrix.iloc[count])
    tax_family = list(family_dummy_matrix.iloc[count])
    tax_genus = list(genus_dummy_matrix.iloc[count])
    tax_species = list(species_dummy_martix.iloc[count])
    taxonomy = tax_phylum + tax_class + tax_order + tax_family + tax_genus + tax_species
    #contig_table['taxonomy_matrix'][count] = taxonomy #",".join(taxonomy)
    taxonomy_matrix_dict[contig] = taxonomy
    if contig_table['num_single_copies'][count] > 0:
        vizbin_x = contig_table['vizbin_x_x'][count]
        vizbin_y = contig_table['vizbin_y_x'][count]
        length = contig_table['length_x'][count]
        cov = contig_table['cov_x'][count]
        gc = contig_table['gc_x'][count]
        coding_density = contig_table['coding_density'][count]
        #Binary list
        label = known_genome
        #features.append([vizbin_x,vizbin_y,cov,coding_density])
        features.append([vizbin_x,vizbin_y,cov] + taxonomy)
        short_features.append([gc,cov]) #Short classifier
        labels.append(label)
        length_weight.append(length)
        #Reassign "unclustered" to known genome assignment
        #contig_table['ML_expanded_clustering'][count] = known_genome

def jackknife_training(features,labels):
    #Function to randomly subsample data into halves (hence 0.5), train
    #ML-classifier and make prediction. Used iteratively in
    #calculate_bootstap_replicates() function (see below)
    train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.50)
    my_classifier = tree.DecisionTreeClassifier()
    my_classifier = my_classifier.fit(train_features,train_labels)
    predictions = my_classifier.predict(test_features)
    return my_classifier

def calculate_bootstap_replicates(feature_array,iterations = 100):
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
    #print top_prediction,confidence_percent
    return top_prediction,confidence_percent

ML_predictions_dict = {}
ML_recruitment_list = []
prediction_accuracy_list = []
#Recruit unclustered seqquences
print("Recruiting unclustered sequences. This could take a while...")
for count,contig in enumerate(contig_table['contig']):
    taxonomy = taxonomy_matrix_dict[contig]
    vizbin_x_x = contig_table.iloc[count]['vizbin_x_x']
    vizbin_y_x = contig_table.iloc[count]['vizbin_y_x']
    cov_x = contig_table.iloc[count]['cov_x']
    composition_feature_list = [vizbin_x_x,vizbin_y_x,cov_x]
    #Concatenate composition and taxonomy features into single list
    single_np_array = np.array([composition_feature_list + taxonomy])
    contig_length = contig_table.iloc[count]['length_x']
    cluster = contig_table.iloc[count][cluster_column_name]
    #If contig is unclustered, make ML prediction
    if cluster == unclustered_name:
        ML_prediction,confidence = calculate_bootstap_replicates(single_np_array,bootstrap_iterations)
        ML_predictions_dict[contig] = ML_prediction,confidence
        #print("ML predictions and jackknife confidence for contig {}: {},{}".format(contig, ML_prediction,confidence))
        #If it the prediction passes confidence cutoff
        if confidence >= confidence_cutoff:
            #Add prediction to ML_recruitment_list
            ML_recruitment_list.append(ML_prediction)
            #For synthetic,pre-defined community (e.g., "MIX_51")
            if ML_prediction == contig_table['known_genome'][count]:
                prediction_accuracy_list.append(1)
            else:
                prediction_accuracy_list.append(0)
        else:
            ML_recruitment_list.append(unclustered_name)
    else:
        ML_recruitment_list.append(cluster)

num_predictions = len(ML_predictions_dict)
num_confident_predictions = len(prediction_accuracy_list)
print("{} of {} predictions were {}% confident".format(num_confident_predictions,num_predictions,confidence_cutoff))

contig_table['ML_expanded_clustering'] = ML_recruitment_list
contig_table.to_csv("1-Dec-16_MIX_51_sk_learn_jackknife.tab",sep="\t",index=False)
