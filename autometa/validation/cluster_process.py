#!/usr/bin/env python2.7

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

# Program to do post processing after binning

from __future__ import division
import argparse
import os
import pandas as pd
from Bio import SeqIO
import subprocess
import math

def run_command(command_string, stdout_path = None):
	# Function that checks if a command ran properly. If it didn't, then print an error message then quit
	print('cluster_process.py, run_command: ' + command_string)
	if stdout_path:
		f = open(stdout_path, 'w')
		exit_code = subprocess.call(command_string, stdout=f, shell=True)
		f.close()
	else:
		exit_code = subprocess.call(command_string, shell=True)

	if exit_code != 0:
		print('cluster_process.py: Error, the command:')
		print(command_string)
		print('failed, with exit code ' + str(exit_code))
		exit(1)

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

parser = argparse.ArgumentParser(description='Script to summarize the assembly characteristics and taxonomy of binned clusters, as well producing individual cluster fasta files')
parser.add_argument('-b', '--bin_table', metavar='<bin.tab>', help='path to the output from either run_autometa.py or ML_recruitment.py', required=True)
parser.add_argument('-c', '--column', metavar='<bin column name>', help='the name of the column to use for binning purposes', default='cluster')
parser.add_argument('-f', '--fasta', metavar='<assembly.fasta>', help='path to the assembly used to make the bin table', required=True)
parser.add_argument('-o', '--output_dir', metavar='<dir>', help='path to the directory where output files will go', default='.')
parser.add_argument('-k', '--kingdom', metavar='<archaea|bacteria>', help='kingdom to consider', choices=['bacteria', 'archaea'], default='bacteria')
parser.add_argument('-t', '--do_taxonomy', help='carry out taxonomic analysis on the clusters (you must have already run make_taxonomy_table.py)', action='store_true')
parser.add_argument('-db', '--db_dir', metavar='<dir>', help='Path to directory with taxdump files')
args = vars(parser.parse_args())

bin_table_path = args['bin_table']
cluster_column_heading = args['column']
fasta_path = args['fasta']
output_dir = args['output_dir']
kingdom = args['kingdom']
db_dir = args['db_dir']
do_taxonomy = args['do_taxonomy']

# Check paths exist
if not os.path.isfile(bin_table_path):
	print('Error! Could not find a bin table at the following path: ' + bin_table_path)
	exit(1)

if not os.path.isfile(fasta_path):
	print('Error! Cannot find a fasta file at the following path: ' + fasta_path)
	exit(1)

# If the user has specified --do_taxonomy, then they also need to specify --db_dir
if do_taxonomy:
	if not db_dir:
		print('Error! If you want to analyze taxonomy, you need to specify a path to database files (--db_dir)')
		exit(1)

	if not os.path.isdir(db_dir):
		print('Error! DB dir ' + db_dir + ' does not exist')
		exit(1)
	if not os.path.isfile(db_dir + '/names.dmp'):
		print('Error! Cannot find names.dmp in ' + db_dir)
		exit(1)

	if not os.path.isfile(db_dir + '/nodes.dmp'):
		print('Error! Cannot find nodes.dmp in ' + db_dir)
		exit(1)

# Make output directory if it isn't already there
if not os.path.isdir(output_dir):
	os.makedirs(output_dir)

master_table = pd.read_table(bin_table_path)

# Format check for the table
columns_to_check = [cluster_column_heading, 'contig', 'length', 'cov', 'single_copy_PFAMs']
if do_taxonomy:
	columns_to_check.append('taxid')

for column in columns_to_check:
	if column not in master_table.columns:
		print('Error! Could not find a column called ' + column + ' in table ' + bin_table_path)
		exit(1)

contig_info = dict() # Will hold dictionaries for each contig, storing length, gc and cov (to calculate weighted av. of gc and cov later for each cluster)
cluster_contigs = dict() # Stores the cluster for each contig
markers_in_cluster = dict() # Keyed by cluster and then PFAM, stores number of the PFAM in each cluster

for i,row in master_table.iterrows():
	contig = row['contig']
	length = int(row['length'])
	cov = float(row['cov'])
	gc = float(row['gc'])
	cluster = row[cluster_column_heading]

	contig_info[contig] = { 'length': length, 'cov': cov, 'gc': gc }
	cluster_contigs[contig] = cluster

	if cluster not in markers_in_cluster:
		markers_in_cluster[cluster] = dict()

	# Protect for instances where single_copy_PFAMs is empty, and gets converted to nan (a float) by pandas
	if isinstance(row['single_copy_PFAMs'], float):
		continue

	pfam_list = row['single_copy_PFAMs'].split(',')

	for pfam in pfam_list:
		if pfam not in markers_in_cluster[cluster]:
			markers_in_cluster[cluster][pfam] = 1
		else:
			markers_in_cluster[cluster][pfam] += 1

# Load fasta file using biopython, and split into clusters
cluster_sequences = dict() # Keyed by cluster, will hold lists of seq objects
for seq_record in SeqIO.parse(fasta_path, 'fasta'):
	seq_name = str(seq_record.id)
	if seq_name in cluster_contigs:
		cluster = cluster_contigs[seq_name]
	else:
		continue

	if cluster not in cluster_sequences:
		cluster_sequences[cluster] = list()
	cluster_sequences[cluster].append(seq_record)

# Output summary table plus individual fasta files
summary_table_path = output_dir + '/cluster_summary_table'
summary_table = open(summary_table_path, 'w')
summary_table.write('cluster\tsize\tlongest_contig\tn50\tnumber_contigs\tcompleteness\tpurity\tav_cov\tav_gc\n')

for cluster in cluster_sequences:
	attributes = assess_assembly(cluster_sequences[cluster])
	total_size = attributes['size']
	longest_contig = attributes['largest_sequence']
	n50 = attributes['n50']
	number_contigs = attributes['number_sequences']

	if kingdom == 'bacteria':
		total_markers = 139
	elif kingdom == 'archaea':
		total_markers = 162

	number_markers_found = 0
	number_unique_markers = len(markers_in_cluster[cluster])
	for pfam in markers_in_cluster[cluster]:
		number_markers_found += markers_in_cluster[cluster][pfam]

	# The following protects for the edge case where a cluster has zero marker genes
	if number_unique_markers == 0:
		completness = 'unknown'
	else:
		completeness = (number_unique_markers / total_markers) * 100
	if number_markers_found == 0:
		purity = 'unknown'
	else:
		purity = (number_unique_markers / number_markers_found) * 100

	# Calculate average GC and cov, weighted by sequence length
	weighted_gc_av = 0.0
	weighted_cov_av = 0.0
	for seq_record in cluster_sequences[cluster]:
		seq_name = str(seq_record.id)
		seq_length = contig_info[seq_name]['length']
		seq_gc = contig_info[seq_name]['gc']
		seq_cov = contig_info[seq_name]['cov']
		seq_length_frac = seq_length / total_size
		weighted_gc_av += seq_gc * seq_length_frac
		weighted_cov_av += seq_cov * seq_length_frac

	# Write line in summary table
	output_string = '\t'.join([cluster, str(total_size), str(longest_contig), str(n50), str(number_contigs), str(completeness), str(purity), str(weighted_cov_av), str(weighted_gc_av)])
	summary_table.write(output_string + '\n')

	# Write individual fasta file
	fasta_output_path = output_dir + '/cluster_' + cluster + '.fasta'
	SeqIO.write(cluster_sequences[cluster], fasta_output_path, 'fasta')

summary_table.close()

if do_taxonomy:
	# Now run cluster_taxonomy.py
	taxonomy_output_path = output_dir + '/cluster_taxonomy.tab'
	cluster_taxonomy_command = 'cluster_taxonomy.py -t {} -c {} -x {} -o {}'.format(bin_table_path, cluster_column_heading, db_dir, taxonomy_output_path)

	run_command(cluster_taxonomy_command)
