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

# Docker wrapper for make_taxonomy_table.py
# Aims to be invisible to the user (use in exactly the same way as make_taxonomy_table.py)

import argparse
import os
import subprocess
import sys

def run_command(command_string, stdout_path = None):
	# Function that checks if a command ran properly. If it didn't, then print an error message then quit
	print('make_taxonomy_table_docker.py, run_command: ' + command_string)
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

pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)

#argument parser
parser = argparse.ArgumentParser(description="Script to generate the contig taxonomy table.", epilog="Output will be directed to recursive_dbscan.py")
parser.add_argument('-a', '--assembly', metavar='<assembly.fasta>', help='Path to metagenomic assembly fasta', required=True)
parser.add_argument('-p', '--processors', metavar='<int>', help='Number of processors to use.', type=int, default=1)
parser.add_argument('-db', '--db_dir', metavar='<dir>', help='Path to directory with taxdump, protein accessions and diamond (NR) protein files. If this path does not exist, will create and download files.', required=False, default=autometa_path + '/databases')
parser.add_argument('-l', '--length_cutoff', metavar='<int>', help='Contig length cutoff to consider for binning in bp', default=10000, type = int)
parser.add_argument('-u', '--update', required=False, action='store_true',\
 help='Checks/Adds/Updates: nodes.dmp, names.dmp, accession2taxid, nr.dmnd files within specified directory.')
parser.add_argument('-v', '--cov_table', metavar='<coverage.tab>', help="Path to coverage table made by calculate_read_coverage.py. If this is not specified then coverage information will be extracted from contig names (SPAdes format)", required=False)
parser.add_argument('-o', '--output_dir', metavar='<dir>', help='Path to directory to store output files', default='.')
parser.add_argument('-s', '--single_genome', help='Specifies single genome mode', action='store_true')

args = vars(parser.parse_args())

db_dir_path = args['db_dir'].rstrip('/')
num_processors = args['processors']
length_cutoff = args['length_cutoff']
fasta_path = args['assembly']
fasta_assembly_prefix = os.path.splitext(os.path.basename(args['assembly']))[0]
prodigal_output = fasta_assembly_prefix + "_filtered.orfs"
prodigal_diamond = prodigal_output + ".blastp"
#add_contig_path = pipeline_path
filtered_assembly = fasta_assembly_prefix + "_filtered.fasta"
cov_table = args['cov_table']
update = args['update']
output_dir = args['output_dir']
single_genome_mode = args['single_genome']

#check if fasta in path
if not os.path.isfile(fasta_path):
	print "Could not find {}...".format(fasta_path)
	logger.debug('Could not find {}...'.format(fasta_path))
	exit(1)

# We have to create the output dir if it doesn't exist
if not os.path.isdir(output_dir):
	os.makedirs(output_dir)

# If db_dir is given and it doesn't exist, then it must be created
if db_dir_path:
	if not os.path.isdir(db_dir_path):
		os.makedirs(db_dir_path)
db_dir_path_absolute = os.path.abspath(db_dir_path)

#If the assembly fasta is not already in the output dir, we need to copy it there so that the docker container can see it
output_dir_absolute = os.path.abspath(output_dir)
fasta_path_absolute = os.path.abspath(fasta_path)
fasta_directory = '/'.join(fasta_path_absolute.split('/')[:-1])
fasta_filename = fasta_path_absolute.split('/')[-1]

if output_dir_absolute != fasta_directory:
	# This means we need to copy the fasta assembly to the output directory
	run_command('cp ' + fasta_path_absolute + ' ' + output_dir_absolute + '/')

#If a coverage table is give, it must already exist.
#If it exists, then it also should be copied to the output directory if it isn't already there
if cov_table:
	if not os.path.isfile(cov_table):
		print('Error! Could not find coverage table at the following path: ' + cov_table)
		exit(1)
	else:
		cov_table_absolute = os.path.abspath(cov_table)
		cov_table_directory = '/'.join(cov_table_absolute.split('/')[:-1])
		cov_table_filename = cov_table_absolute.split('/')[1]
		if output_dir_absolute != cov_table_directory:
			run_command('cp ' + cov_table_absolute + ' ' + output_dir_absolute + '/')

# Construct make_taxonomy_table.py command to pass to the docker container
make_taxonomy_table_command = 'make_taxonomy_table.py --assembly /output/{} --processors {} --db_dir /databases --length_cutoff {} --output_dir /output'.format(\
	fasta_filename, num_processors, length_cutoff)

if cov_table:
	make_taxonomy_table_command = make_taxonomy_table_command + ' --cov_table {}'.format(cov_table_filename)

if update:
	make_taxonomy_table_command = make_taxonomy_table_command + ' --update'

if single_genome_mode:
	make_taxonomy_table_command = make_taxonomy_table_command + ' --single_genome'
	
# Construct Docker run command
docker_command = 'docker run --volume {}:/output:rw --volume {}:/databases:rw --detach=false --rm jasonkwan/autometa:latest'.format(output_dir_absolute, db_dir_path_absolute)

# Incorporate autometa command
docker_command = docker_command + ' ' + make_taxonomy_table_command

# Run docker container
run_command(docker_command)
