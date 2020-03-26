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

# Docker wrapper for cluster_process.py
# Aims to be invisible to the user (use in exactly the same way as cluster_process.py)

import argparse
import os
import subprocess

def run_command(command_string, stdout_path = None):
	# Function that checks if a command ran properly. If it didn't, then print an error message then quit
	print('cluster_process_docker.py, run_command: ' + command_string)
	if stdout_path:
		f = open(stdout_path, 'w')
		exit_code = subprocess.call(command_string, stdout=f, shell=True)
		f.close()
	else:
		exit_code = subprocess.call(command_string, shell=True)

	if exit_code != 0:
		print('run_autometa_docker.py: Error, the command:')
		print(command_string)
		print('failed, with exit code ' + str(exit_code))
		exit(1)

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
	else:
		if not os.path.isfile(db_dir + '/names.dmp'):
			print('Error! Cannot find names.dmp in ' + db_dir)
			exit(1)

		if not os.path.isfile(db_dir + '/nodes.dmp'):
			print('Error! Cannot find nodes.dmp in ' + db_dir)
			exit(1)

# Make output directory if it isn't already there
if not os.path.isdir(output_dir):
	os.makedirs(output_dir)

output_dir_absolute = os.path.abspath(output_dir)
db_dir_absolute = os.path.abspath(db_dir)

# If fasta and bin table are not already in the output dir, we should copy them there now
fasta_path_absolute = os.path.abspath(fasta_path)
fasta_dir = '/'.join(fasta_path_absolute.split('/')[:-1])
fasta_filename = fasta_path_absolute.split('/')[-1]

if fasta_dir != output_dir_absolute:
	run_command('cp ' + fasta_path_absolute + ' ' + output_dir_absolute + '/')

bin_table_path_absolute = os.path.abspath(bin_table_path)
bin_table_dir = '/'.join(bin_table_path_absolute.split('/')[:-1])
bin_table_filename = bin_table_path_absolute.split('/')[-1]

if bin_table_dir != output_dir_absolute:
	run_command('cp ' + bin_table_path_absolute + ' ' + output_dir_absolute + '/')

# Construct cluster_process.py command
cluster_process_command = 'cluster_process.py --bin_table /output/{} --column {} --fasta /output/{} --output_dir /output --kingdom {}'.format(\
	bin_table_filename, cluster_column_heading, fasta_filename, kingdom)

if do_taxonomy:
	cluster_process_command = cluster_process_command + ' --do_taxonomy --db_dir /databases'

# Construct docker run command
docker_command = 'docker run --volume {}:/output:rw --detach=false --rm'.format(output_dir_absolute)

if db_dir:
	docker_command = docker_command + ' --volume {}:/databases:rw'.format(db_dir_absolute)

docker_command = docker_command + ' jasonkwan/autometa:latest'

# Incorporate autometa command
docker_command = docker_command + ' ' + cluster_process_command

# Run docker container
run_command(docker_command)
