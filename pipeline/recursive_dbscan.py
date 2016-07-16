#!/usr/bin/env python

# Program to carry out secondary dbscan clustering on a re-run of vizbin on unclustered contigs (using low "perplexity" setting)
# When you run vizbin this way, Davies-Bouldin index seems to fail because the groups are quite diffuse, although DBSCAN seems to do
# fine with the correct eps value.  Here we judge DBSCAN results based on the number of pure clusters (sum of total purity)

import rpy2.robjects as robjects # Bridge to R code
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
import sys
import os
import pandas as pd
import csv
import argparse
import copy
import numpy
from Bio import SeqIO
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

pandas2ri.activate()
pp = pprint.PrettyPrinter(indent=4)

# Load R libaries
rbase = importr('base')
dbscan = importr('dbscan')

# Define R functions
robjects.r('''
	get_table <- function(path) {
		input_data <- read.table(path, header=TRUE)
		return (input_data)
	}
	dbscan_simple <- function(input_data_frame, eps) {
		# Funny things happen if there is already a 'db_cluster' column, let's delete it!
		if ("db_cluster" %in% colnames(input_data_frame))
		{
			input_data_frame$db_cluster <- NULL
		}

		d <- data.frame(input_data_frame$vizbin_x, input_data_frame$vizbin_y)

		db <- dbscan(d, eps=eps, minPts=3)
		output_table <- data.frame(input_data_frame, db_cluster = db$cluster )

		return(output_table)
	}
''')

dbscan_simple = robjects.r['dbscan_simple']
get_table = robjects.r['get_table']

parser = argparse.ArgumentParser(description="Prototype script to automatically carry out secondary clustering of vizbin coordinates based on DBSCAN and cluster purity")
parser.add_argument('-m','--marker_tab', help='Output of make_marker_table.py', required=True)
#parser.add_argument('-v','--vizbin_tab', help='Table containing vizbin coordinates', required=True)
parser.add_argument('-d','--domain', help='Microbial domain (bacteria|archaea)', default='bacteria')
parser.add_argument('-f','--fasta', help='Assembly FASTA file', required=True) 
parser.add_argument('-o','--outdir', help='Path of directory for output', required=True)
parser.add_argument('-c','--contigtable', help='Output of make_contig_table.py', required=True)
parser.add_argument('-p','--processors', help='Number of processors used', default=1)
#parser.add_argument('-p','--purity_cutoff', help='Cutoff (%) used to count number of pure clusters', default=90)
#parser.add_argument('-c','--completeness_cutoff', help='Cutoff (%) used to count number of complete clusters', default=90)
args = vars(parser.parse_args())

hmm_table_path = args['marker_tab']
#vizbin_table_path = args['vizbin_tab']
domain = args['domain']
fasta_path = args['fasta']
outdir = os.path.abspath(args['outdir'])
contig_table = args['contigtable']
processors = args['processors']
#purity_cutoff = int(args['purity_cutoff'])
#completeness_cutoff = int(args['completeness_cutoff'])

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
		total_unique = 0
		for marker in marker_totals[cluster]:
			if marker_totals[cluster][marker] == 1:
				total_unique += 1
		completeness = (float(total_markers) / expected_number) * 100
		# Protect from divide by zero
		if total_markers == 0:
			purity = 0
		else:
			purity = (float(total_unique) / total_markers) * 100
		cluster_details[cluster] = { 'completeness': completeness, 'purity': purity }

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

def runDBSCANs(r_table):
	# Carry out DBSCAN, starting at eps=0.3 and continuing until there is just one group
	current_eps = 0.3
	db_tables = {} # Will be keyed by eps
	number_of_clusters = float('inf')
	while(number_of_clusters > 1):
		#print ('current eps: ' + str(current_eps))
		dbscan_output_r = dbscan_simple(r_table, current_eps)
		dbscan_output_pd = pandas2ri.ri2py(dbscan_output_r)
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

	# For impure clusters, output vizbin table
	# First, find pure clusters
	best_db_table = table_dictionary[best_eps_value]

	complete_and_pure_clusters = {}
	other_clusters = {}
	for cluster in cluster_info[best_eps_value]:
		completeness = cluster_info[best_eps_value][cluster]['completeness']
		purity = cluster_info[best_eps_value][cluster]['purity']

		# Explicitly remove the noise cluster (0)
		if cluster == 0:
			other_clusters[cluster] = 1
			continue

		if completeness > completeness_cutoff and purity > purity_cutoff:
			complete_and_pure_clusters[cluster] = 1
		else:
			other_clusters[cluster] = 1

	# Subset the data frame
	subset_other_db_table = best_db_table
	for cluster in complete_and_pure_clusters:
		subset_other_db_table = subset_other_db_table[subset_other_db_table['db_cluster'] != cluster]

	# We now make an r table from subset_other_db_table, with the db_cluster column stripped
	subset_other_db_table = subset_other_db_table.drop('db_cluster', 1)
	unclustered_r = pandas2ri.py2ri(subset_other_db_table)

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

	return output_cluster_info, output_contig_cluster, unclustered_r

def run_VizBin(fasta, output_filename):
	if os.path.isfile(output_filename):
		print "VizBin output already exists!"
		print "Continuing to next step..."
		logger.info('VizBin output already exists!')
		logger.info('Continuing to next step...')
		return None
	else:
		print "Runnign k-mer based binning..."
		logger.info('Running k-mer based binning...')
		subprocess.call("java -jar {}VizBin-dist.jar -i {} -o points.txt -t {}".format(autometa_path + "/VizBin/dist/", fasta, processors), shell = True, stderr=subprocess.STDOUT)
		#tmp_path = subprocess.check_output("ls /tmp/map* -dlt | grep {} | head -n1".format(username), shell=True).rstrip("\n").split()[-1]
		return fasta

def process_and_clean_VizBin(input_fasta,contig_table,output_table_name):	
	subprocess.call("{}/vizbin_process.pl {} points.txt > vizbin_table.txt".format(pipeline_path, input_fasta), shell=True) # Note: we assume here that the cutoff is above 100 bp
	subprocess.call("tail -n +2 vizbin_table.txt | sort -k 1,1 > vizbin_table_sort.txt", shell=True)
	subprocess.call("tail -n +2 {} | sort -k 1,1 > contig_table_sort.txt".format(contig_table), shell=True)
	subprocess.call("head -n 1 vizbin_table.txt > vizbin_header.txt", shell=True)
	subprocess.call("head -n 1 {} > contig_header.txt".format(contig_table), shell=True)
	subprocess.call("join contig_header.txt vizbin_header.txt | sed $'s/ /\t/g' > joined_header.txt", shell=True)
	subprocess.call("cat joined_header.txt > {}".format(output_table_name), shell=True)
	subprocess.call("join contig_table_sort.txt vizbin_table_sort.txt | sed $'s/ /\t/g' >> {}".format(output_table_name), shell=True)
		#subprocess.call("touch {}; cat joined_header.txt >> {}; join contig_table_sort.txt vizbin_table_sort.txt | cat >> {}; sed $'s/ /\t/g' contig_vizbin.tab", shell=True,stdout=FNULL, stderr=subprocess.STDOUT)
		#Delete most recent /tmp/map* directory if it's the same user 
		#subprocess.call("rm -rf {}".format(tmp_path), shell = True)
		#need more elegant way to clean up
	subprocess.call("rm *.txt", shell = True)



start_time = time.time()
FNULL = open(os.devnull, 'w')

pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)
username = getpass.getuser()



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

# Make a master table that will be updated as clustering is done




#pp.pprint(db_tables)
#pdb.set_trace()

# For the first few rounds, we are going to use completeness_cutoff = 90 and purity_cutoff = 90
global_cluster_info = {}
global_cluster_contigs = {}
#completeness_cutoff = 90
completeness_cutoff = 20
purity_cutoff = 90
round_counter = 0
current_fasta = fasta_path
vizbin_counter = 0
master_table = None

# Load fasta sequences
assembly_seqs = {}
for seq_record in SeqIO.parse(fasta_path, 'fasta'):
	assembly_seqs[seq_record.id] = seq_record

while True:
	# Run vizbin
	vizbin_counter += 1
	logger.info('Running VizBin round ' + str(vizbin_counter))
	current_vizbin_output = 'vizbin' + str(vizbin_counter) + '.tab'
	# Carry out first vizbin run
	process_and_clean_VizBin(run_VizBin(current_fasta,current_vizbin_output),contig_table, current_vizbin_output)

	abs_vizbin_path = os.path.abspath(current_vizbin_output)
	vizbin_r = get_table(abs_vizbin_path)

	if vizbin_counter == 1:
		master_table = pandas2ri.ri2py(vizbin_r)

	# Output current vizbin table
	vizbin_output_path = 'vizbin' + str(vizbin_counter) + '.tab'
	vizbin_pd = pandas2ri.ri2py(vizbin_r)
	vizbin_pd.to_csv(path_or_buf=vizbin_output_path, sep="\t", index=False, quoting=csv.QUOTE_NONE)

	local_vizbin_round = 0
	while True:
		round_counter += 1
		local_vizbin_round += 1
		logger.info('Running DBSCAN round ' + str(round_counter))
		db_tables = runDBSCANs(current_r_table)
		cluster_information, contig_cluster_dictionary, unclustered_r_table = assessDBSCAN(db_tables, contig_markers, domain, completeness_cutoff, purity_cutoff)
		current_r_table = unclustered_r_table

		if not cluster_information:
			break

		# Populate the global data structures
		for	cluster in cluster_information:
			new_cluster_name = 'vizbin' + str(vizbin_counter) + '_round' + str(round_counter) + '_' + str(cluster)
			global_cluster_info[new_cluster_name] = cluster_information[cluster]

		for contig in contig_cluster_dictionary:
			new_cluster_name = 'vizbin' + str(vizbin_counter) + '_round' + str(round_counter) + '_' + str(contig_cluster_dictionary[contig])
			global_cluster_contigs[contig] = new_cluster_name

	if local_vizbin_round == 1:
		# This means that after the last vizbin run, dbscan only ran once, meaning that no clusters were found upon first run, and we are done
		break

	# If we are not done, write a new fasta for the next vizbin run
	# First convert the r table to a pandas table
	unclustered_pd = pandas2ri.ri2py(current_r_table)
	unclustered_seqrecords = []
	for i, row in unclustered_pd.iterrows():
		contig = row['contig']
		unclustered_seqrecords.append(assembly_seqs[contig])

	current_fasta = 'unclustered_for_vizbin' + str(vizbin_counter + 1) + '.fasta'
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









