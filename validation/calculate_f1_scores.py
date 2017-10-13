#!/usr/bin/env python

#1. Compile data structure, for contig each contig, record: length, cluster, reference genome
#2. Compile data structure for each reference genome (list of contigs, list of clusters?)

import pandas as pd
import operator
import argparse

parser = argparse.ArgumentParser(description='Script to Calculate F1 from master\
 contig table (for heat map).')
parser.add_argument('-c','--cluster_column', help='Name of the cluster column', required=True)
parser.add_argument('-t','--contig_tab', help='Name of master contig table file', required=True)
parser.add_argument('-u','--unclustered_term', help='Name of unclustered bin. Default = "unclustered"', default="unclustered")
parser.add_argument('-o','--out', help='Tab delimited output table with F1 scores', default="F1.tab")
args = vars(parser.parse_args())

#args = {"cluster_column":"cluster","contig_tab":"MIX_51_master_w_taxonomy.tab","unclustered_term":0, "out":"MIX_51_F1_dbscan_test.txt"}

contig_table = pd.read_csv(args['contig_tab'], sep = "\t")
contig_table[args['cluster_column']].astype(basestring)

#1. Compile data structures, for each reference genome and each contig (length, cluster)
reference_genome_dict = {}
contig_dict = {}
global_cluster_dict = {}
for index,row in contig_table.iterrows():
    contig = row['contig']
    known_genome = row['reference_genome']
    length = row['length']
    cluster = row[args['cluster_column']]
    #build reference genome dict
    if known_genome not in reference_genome_dict:
        reference_genome_dict[known_genome] = [contig]
    else:
        reference_genome_dict[known_genome].append(contig)
    #build contig dict
    contig_dict[contig] = {}
    contig_dict[contig]['length'] = length
    contig_dict[contig]['cluster'] = cluster
    #build global cluster dict
    if isinstance(cluster,float):
        cluster = str(cluster).split(".")[0]
    if str(cluster) not in global_cluster_dict:
        global_cluster_dict[str(cluster)] = [contig]
    else:
        global_cluster_dict[str(cluster)].append(contig)

#2. Calculate precision and recall

def find_dominant_cluster_length(contig_list,contig_dict,reference_genome):
    #Function to calculate recall based on contig_list and contig_dict
    #A. Identify which cluster best/most represents reference genome
    cluster_dict = {}
    for contig in contig_list:
        cluster = contig_dict[contig]['cluster']
        length = contig_dict[contig]['length']
        if isinstance(cluster,float):
            cluster = str(cluster).split(".")[0]
        if str(cluster) not in cluster_dict:
            cluster_dict[str(cluster)] = length
        else:
            cluster_dict[str(cluster)] += length

    #Find dominant cluster
    max_cluster_length = max(cluster_dict.values())
    #make sure that there aren't two clusters with the same length
    try:
        assert cluster_dict.values().count(max_cluster_length) <= 1
    except AssertionError:
        num_clusters_with_same_mapped_length = cluster_dict.values().count(max_cluster_length)
        print
        #I don't think this would actually affect the F1 calculation, either case it would be the same
        print("Warning: There are {} clusters that share a max aligned length of {} for {}".format(num_clusters_with_same_mapped_length,max_cluster_length,reference_genome))
        dominant_clusters = sorted(cluster_dict.items(), key=operator.itemgetter(1))[-num_clusters_with_same_mapped_length:]
        #print cluster_dict
        print dominant_clusters
        print
    #convert nan to string
    dominant_cluster = str(max(cluster_dict.iteritems(), key=operator.itemgetter(1))[0])
    #if the "dominant" cluster is unclustered, then use next in line, if available
    if dominant_cluster == args['unclustered_term'] and len(cluster_dict.keys()) > 1:
        #Sort into a list of tuple, ascending based on values in dict, grab the second of last, and report the key
        dominant_cluster = sorted(cluster_dict.items(), key=operator.itemgetter(1))[-2][0]
        print("Dominant cluster in '{}', using second largest cluster ({}) instead.".format(args['unclustered_term'],dominant_cluster))

    #return dominant_cluster and length of contigs in that cluster that have that ref genome
    #print dominant_cluster
    #print cluster_dict
    return str(dominant_cluster),cluster_dict[dominant_cluster]

#single test case
#find_dominant_cluster_length(reference_genome_dict['GCA_000466525'],contig_dict)
#print global_cluster_dict
#Now for all the reference genomes
with open(args['out'],"w") as outfile:
    outfile.write("ref_genome\tref_genome_len\tdominant_cluster\ttotal_dominant_cluster_len\tref_aln_dominant_cluster_len\trecall\tprecision\tF1\n")
    print("ref_genome\tref_genome_len\tdominant_cluster\ttotal_dominant_cluster_len\tref_aln_dominant_cluster_len\trecall\tprecision\tF1\n")
    for reference_genome,contig_list in reference_genome_dict.items():
        dominant_cluster,dominant_cluster_length = find_dominant_cluster_length(contig_list,contig_dict,reference_genome)
        reference_genome_length = sum([contig_dict[contig]['length'] for contig in contig_list])
        if dominant_cluster == args['unclustered_term']:
            recall = 0.0
            precision = 0.0
            total_cluster_len = 0.0
            dominant_cluster_length = 0.0
            F1 = 0.0
        else:
            #recall is % of which that dominant cluster recovers the reference genome
            recall = round(dominant_cluster_length/float(reference_genome_length)*100,5)
            #precision is the  % (in length) of dominant cluster that's represented by corresponding reference genome
            total_cluster_len = 0
            ref_genome_len = 0
            for contig in global_cluster_dict[dominant_cluster]:
                contig_length = contig_dict[contig]['length']
                #If it's in the list of ref genome contigs, add to precision count
                if contig in contig_list:
                    ref_genome_len += contig_length
                #else: print("\t{} is in bin: {}, but does not have ref genome: {}".format(contig,dominant_cluster,reference_genome))
                total_cluster_len += contig_length
            precision = round(ref_genome_len/float(total_cluster_len)*100)
            F1 = round(2 * (precision*recall)/(precision+recall),5)

        print reference_genome,reference_genome_length,dominant_cluster,total_cluster_len,dominant_cluster_length,recall,precision,F1
        out_list = [reference_genome,reference_genome_length,dominant_cluster,total_cluster_len,dominant_cluster_length,recall,precision,F1]
        outfile.write("\t".join(map(str,out_list)) + "\n")
    outfile.write("\n")
