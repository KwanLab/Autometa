#!/usr/bin/env python

# Program that determines the completeness of clusters called by dbscan

import sys
import argparse 
from Bio import SeqIO
import pdb
import pandas as pd
import numpy as np

def assess_assembly(seq_record_list):
	assembly_size = sum(len(seq) for seq in seq_record_list)
	number_of_sequences = len(seq_record_list)
	sorted_seqs = sorted(seq_record_list, key=len)
	largest_sequence_length = len(sorted_seqs[-1])
	sequence_total = 0
	n50 = None
	for i,seq_object in enumerate(sorted_seqs):
		sequence_total += len(seq_object)
		if sequence_total > (float(assembly_size)/2):
			n50 = len(seq_object)
			break
	return { 'size': assembly_size, 'number_sequences': number_of_sequences, 'largest_sequence': largest_sequence_length, 'n50': n50 }

def pick_low_hanging_fruit(summary_table, marker_table_rows):
	#First, parse table
	summary_table = open(summary_table_path, 'r')
	summary_table_rows = summary_table.read().splitlines()
	summary_table.close()
	
	best_clusters = []
	best_purity = []
	best_completeness = []
	#Go through each cluster and check for requried completeness and requried purity
	for i,line in enumerate(summary_table_rows):
	       if i > 0:
	            line_list = line.split('\t')
	            cluster = line_list[0]
	            completeness = line_list[5]
	            purity = line_list[6]

	             #Checking each cluster
	            if (float(completeness) >= float(completeness_min)) and (float(purity) >= float(purity_min)):
	                    best_clusters.append(cluster)
	                    best_completeness.append(completeness)
	                    best_purity.append(purity)
	#Printing out the clusters that passed
	#print "\n%d out of %d clusters are considered to be 'low hanging fruit':" % (len(best_clusters), (len(summary_table_rows)-1))
	#for i, item in enumerate(best_clusters):
	#	print "Cluster: " + best_clusters[i] + " with completeness: " + best_completeness[i] + " and purity: " + best_purity[i]
	#print "\n"
	return best_clusters

parser = argparse.ArgumentParser(description='Script to determine the completeness of clusters called by dbscan')
parser.add_argument('-d','--dbscantable', help='table containing dbscan information', required=True)
parser.add_argument('-c','--column', help='bin column name in dbscan table', default = 'db.cluster')
parser.add_argument('-m','--markertable', help='marker table created with make_marker_table', required=True)
parser.add_argument('-f','--fasta', help='contig fasta file', required=True)
parser.add_argument('-o','--output', help='output directory for summary table and cluster fasta files', required=True)
parser.add_argument('-k','--kingdom', help='kingdom (bacteria|archaea)', default = 'bacteria')
parser.add_argument('-cc','--cluster_completeness', help='fasta files outputed determined by cluster_completeness', default = 20)
parser.add_argument('-cm', '--completeness_minimum',help="Minimum completeness required to be considered 'low hanging fruit'", default=80)
parser.add_argument('-p', '--purity_minimum',help="Minimum purity required to be considered 'low hanging fruit'", default=90)
args = vars(parser.parse_args())

dbscan_table_path = args['dbscantable']
completeness_min = args['completeness_minimum']
purity_min = args['purity_minimum']
cluster_column_heading = args['column']
marker_table_path = args['markertable']
fasta_file_path = args['fasta']
output_prefix = args['output']
kingdom = args['kingdom']
cluster_completeness = float(args['cluster_completeness'])
# Input varification *TO DO*
# Check paths exist
# Check that kingdom is either 'bacteria' or 'archaea'

# First go through marker table
marker_table = open(marker_table_path, 'r')
marker_table_rows = marker_table.read().splitlines()
marker_table.close
contig_markers = {}
markers_in_contig = {} # Keyed by contigs, holds total of each marker found
for i,line in enumerate(marker_table_rows):
	if i > 0:
		line_list = line.split('\t')
		pfam_list = line_list[1].split(',')
		contig_name = line_list[0]
		num_single_copies = line_list[2]
<<<<<<< HEAD
		markers_in_contig[contig_name] = num_single_copies
=======
		contig[contig_name] = num_single_copies
>>>>>>> 9d4a1a1a883ddfbd8b2234b9ba1ad3024ac58d96
		if contig_name not in contig_markers:
			contig_markers[contig_name] = {}

		for pfam in pfam_list:
			if pfam in contig_markers[contig_name]:
				contig_markers[contig_name][pfam] += 1
			else:
				contig_markers[contig_name][pfam] = 1


# Now go through dbscan table
dbscan_table = open(dbscan_table_path, 'r')
dbscan_table_rows = dbscan_table.read().splitlines()
dbscan_table.close

dbscan_header_line = dbscan_table_rows[0]
dbscan_header_list = dbscan_header_line.split('\t')
cluster_index = None
contig_index = None
cov_index = None
gc_index = None
cluster_column_found = 0
contig_column_found = 0
for i, heading in enumerate(dbscan_header_list):
	if heading == cluster_column_heading:
		cluster_index = i
		cluster_column_found += 1
	if heading == 'contig':
		contig_index = i
		contig_column_found += 1
	#New
	if heading == 'gc':
		gc_index = i
	if heading == 'cov':
		cov_index = i

if cluster_index is None:
	print 'Error, could not find column ' + cluster_column_heading + ' in dbscan table ' + dbscan_table_path
	sys.exit(2)

if cluster_column_found > 1:
	print 'Error, multiple columns called ' + cluster_column_heading + ' found in ' + dbscan_table_path
	sys.exit(2)

if contig_index is None:
	print 'Error, could not find contig column in ' + dbscan_table_path
	sys.exit(2)

if contig_column_found > 1:
	print 'Error, multiple contig columns found in ' + dbscan_table_path
	sys.exit(2)


markers_in_cluster = {} #Keyed by cluster, stores number of markers
cluster_contigs = {} # Keyed by contig, stores the cluster of each contig
gc_in_cluster = {} # Keyed by cluster, holds total gc
cov_in_cluster = {}  # keyed by cluster, holds total cov
count = {}
stdev_gc = {}
stdev_cov = {}

for i,line in enumerate(dbscan_table_rows):
	if i > 0:
		line_list = line.split('\t')
		contig = line_list[contig_index]
		gc = float(line_list[gc_index])
		cov = float(line_list[cov_index])
		cluster = line_list[cluster_index]

		if cluster not in gc_in_cluster:
			gc_in_cluster[cluster] = 0
		if cluster not in cov_in_cluster:
			cov_in_cluster[cluster] = 0 
		if cluster not in stdev_cov:
			stdev_cov[cluster] = []
		if cluster not in stdev_gc:
			stdev_gc[cluster] = []

		gc_in_cluster[cluster] += gc  
		cov_in_cluster[cluster] += cov
		stdev_cov[cluster].append(cov)
		stdev_gc[cluster].append(gc)
		cluster_contigs[contig] = cluster

		if cluster not in markers_in_cluster:
			markers_in_cluster[cluster] = {}
		if cluster not in count:
			count[cluster] = 0
		count[cluster] += 1

		#pdb.set_trace()
		if contig != "contig":
			for pfam in contig_markers[contig]:
				if pfam in markers_in_cluster[cluster]:
					markers_in_cluster[cluster][pfam] += contig_markers[contig][pfam]
				else:
					markers_in_cluster[cluster][pfam] = contig_markers[contig][pfam]


# Load fasta file using biopython
# Split into clusters
cluster_sequences = {} # Keyed by cluster, will hold lists of seq objects
for seq_record in SeqIO.parse(fasta_file_path, 'fasta'):
	seq_name = seq_record.id
	seq = seq_record.seq
	cluster = None
	check_markers = int(markers_in_contig[seq_name])
	if seq_name in cluster_contigs:
		cluster = cluster_contigs[seq_name]
	else:
		cluster = 'unclaimed'

	if cluster not in cluster_sequences:
		cluster_sequences[cluster] = []
	cluster_sequences[cluster].append(seq_record)
<<<<<<< HEAD
=======

# Need to total up markers in the 'unclaimed' bin
if 'unclaimed' in cluster_sequences:
	for seq_record in cluster_sequences['unclaimed']:
		if 'unclaimed' not in markers_in_cluster:
			markers_in_cluster['unclaimed'] = {}
		contig = seq_record.id
		for pfam in contig_markers[contig]:
			if pfam in markers_in_cluster['unclaimed']:
				markers_in_cluster['unclaimed'][pfam] += contig_markers[contig][pfam]
			else:
				markers_in_cluster['unclaimed'][pfam] = contig_markers[contig][pfam]

>>>>>>> 9d4a1a1a883ddfbd8b2234b9ba1ad3024ac58d96
# Now go through cluster and output table
summary_table_path = output_prefix + '/summary_table'
summary_table = open(summary_table_path, 'w')
summary_table.write('cluster\tsize\tlongest_contig\tn50\tnumber_contigs\tcompleteness\tpurity\tcov\tstdev_cov\tgc_percent\tstdev_gc\tcompleteness_over_{}\n'.format(cluster_completeness))

for cluster in cluster_sequences:
	if cluster == 'unclaimed':
		continue
		
	attributes = assess_assembly(cluster_sequences[cluster])
	if kingdom == 'bacteria':
		total_markers = 139
	else:
		total_markers = 162

	number_unique_markers = 0
	number_of_markers_found = len(markers_in_cluster[cluster])
	for pfam in markers_in_cluster[cluster]:
		if markers_in_cluster[cluster][pfam] == 1:
			number_unique_markers += 1

	average_gc = float(gc_in_cluster[cluster]/count[cluster])
	average_cov = float(cov_in_cluster[cluster]/count[cluster])
	std_gc = np.std(stdev_gc[cluster], ddof=1)
	std_cov = np.std(stdev_cov[cluster], ddof=1)
	completeness = (float(number_of_markers_found)/total_markers) * 100
	purity = (float(number_unique_markers)/number_of_markers_found) * 100
	
	if completeness >= cluster_completeness:
		summary_table.write(str(cluster) + '\t' + str(attributes['size']) + '\t' + str(attributes['largest_sequence']) + '\t' + str(attributes['n50']) + '\t' + str(attributes['number_sequences']) + '\t{0:.6f}'.format(completeness) + 
		'\t{0:.6f}\t'.format(purity) + str(average_cov) + '\t' + str(std_cov) + '\t' + str(average_gc) + '\t' + str(std_gc) + '\t' + 'True'  + '\n')
       		# Now output the fasta file
       		fasta_output_path = output_prefix + '/cluster_' + cluster + '.fasta'
		SeqIO.write(cluster_sequences[cluster], fasta_output_path, 'fasta')
	else:
		summary_table.write(str(cluster) + '\t' + str(attributes['size']) + '\t' + str(attributes['largest_sequence']) + '\t' + str(attributes['n50']) + '\t' + str(attributes['number_sequences']) +
		'\t{0:.6f}'.format(completeness) + '\t{0:.6f}\t'.format(purity) + str(average_cov) + '\t' + str(std_cov) + '\t' + str(average_gc) + '\t' + str(std_gc) + '\t' + 'False'  + '\n')

#Pick out clusters that have markers that are 80% complete and 90% pure. 
best_clusters = pick_low_hanging_fruit(summary_table,marker_table_rows)

fasta_with_markers = []
fasta_without_markers = []
# Load fasta file using biopython
# Split into two fasta files
for seq_record in SeqIO.parse(fasta_file_path, 'fasta'):
        seq_name = seq_record.id
        seq = seq_record.seq
        cluster_contig = cluster_contigs[seq_name]
        check_markers = int(markers_in_contig[seq_name])
        if check_markers >= 1 and cluster_contig not in best_clusters:
                fasta_with_markers.append(seq_record)
        elif check_markers < 1 and cluster_contig not in best_clusters:
                fasta_without_markers.append(seq_record)

fasta_with_markers_output = open("contigs_with_markers.fasta", "w")
SeqIO.write(fasta_with_markers, fasta_with_markers_output, "fasta")
fasta_with_markers_output.close()

fasta_without_markers_output = open("contigs_without_markers.fasta", "w")
SeqIO.write(fasta_without_markers, fasta_without_markers_output, "fasta")
fasta_without_markers_output.close()
