#!/usr/bin/env python

from __future__ import division
import argparse, collections, pandas as pd

parser = argparse.ArgumentParser(description="Program to determine completeness\
and contamination from the master contig table cluster column.")
parser.add_argument('-t','--contig_tab', help='Master contig table', required=True)
parser.add_argument('-c','--cluster_column', help='Cluster column name', \
default="sklearn_recursive_dbscan")
args = vars(parser.parse_args())

contig_table = pd.read_csv(args['contig_tab'], sep = "\t")
cluster_column = args['cluster_column']

#contig_table = pd.read_csv("ATCC-MSA-1001_master_contig.tab", sep = "\t")
#cluster_column = "sklearn_recursive_dbscan"

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

cluster_dict = calculateClusterStats(contig_table,cluster_column)

print("cluster\tcompleteness\tpurity")
for cluster,info_dictionary in cluster_dict.items():
    print("{}\t{}\t{}").format(cluster,info_dictionary['completeness'],info_dictionary['purity'])
