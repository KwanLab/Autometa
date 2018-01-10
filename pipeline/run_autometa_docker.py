#!/usr/bin/env python

# Docker wrapper for run_autometa.py
# Aims to be invisible to the user (use in exactly the same way as run_autometa.py)

import argparse
import os
import subprocess

def run_command(command_string, stdout_path = None):
	# Function that checks if a command ran properly. If it didn't, then print an error message then quit
	print('calculate_read_coverage.py, run_command: ' + command_string)
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

#argument parser
parser = argparse.ArgumentParser(description="Script to run the Autometa pipeline.",\
 epilog="Please do not forget to cite us. Thank you for using Autometa!",\
  formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-a', '--assembly', metavar='<assembly.fasta>', help='Path to metagenomic assembly fasta', required=True)
parser.add_argument('-p', '--processors', metavar='<int>', help='Number of processors to use', type=int, default=1)
parser.add_argument('-l', '--length_cutoff', metavar='<int>', help='Contig length cutoff to consider for binning in bp', default=10000, type=int)
parser.add_argument('-c', '--completeness_cutoff', metavar='<float>', help='Completeness cutoff (in percent) to use for accepting clusters', type=float, default=20.0)
parser.add_argument('-k', '--kingdom', metavar='<archaea|bacteria>', help='Kingdom to consider',\
choices=['bacteria','archaea'], default = 'bacteria')
parser.add_argument('-t', '--taxonomy_table', metavar='<taxonomy.tab>', help='Path to output of make_taxonomy_table.py')
parser.add_argument('-o', '--output_dir', metavar='<dir>', help='Path to directory to store all output files', default = '.')
parser.add_argument('-r', '--ML_recruitment', help='Use ML to further recruit unclassified contigs', action='store_true')
parser.add_argument('-m', '--maketaxtable', action='store_true',\
help='runs make_taxonomy_table.py before performing autometa binning. Must specify databases directory (-db)')
parser.add_argument('-db', '--db_dir', metavar='<dir>', help="Path to directory with taxdump files. If this doesn't exist, the files will be automatically downloaded", required=False)
parser.add_argument('-v', '--cov_table', metavar='<coverage.tab>', help="Path to coverage table made by calculate_read_coverage.py. If this is not specified then coverage information will be extracted from contig names (SPAdes format)", required=False)

args = vars(parser.parse_args())

length_cutoff = args['length_cutoff']
fasta_assembly = os.path.abspath(args['assembly'])
processors = args['processors']
cluster_completeness = args['completeness_cutoff']
kingdom = args['kingdom'].lower()
taxonomy_table_path = args['taxonomy_table']
output_dir = args['output_dir']
do_ML_recruitment = args['ML_recruitment']
make_tax_table = args['maketaxtable']
db_dir_path = args['db_dir']
cov_table = args['cov_table']

# Argument error checks

#check if appropriate databases specified for make taxonomy table
if make_tax_table and not db_dir_path:
	print("Must specify databases directory (-db)")
	exit(1)

#check if fasta in path
if not os.path.isfile(fasta_assembly):
	print "Could not find {}...".format(fasta_assembly)
	logger.debug('Could not find {}...'.format(fasta_assembly))
	exit(1)

# We have to create the output dir if it doesn't exist
if not os.path.isdir(output_dir):
	os.makedirs(output_dir)

# If db_dir is given and it doesn't exist, then it must be created
if db_dir_path:
	if not os.path.isdir(db_dir_path):
		os.makedirs(db_dir_path)
	db_dir_path_absolute = os.path.abspath(db_dir_path)

#if make_tax_table specified but taxonomy_table_path not defined
if make_tax_table and not taxonomy_table_path:
	taxonomy_table_path = output_dir + '/taxonomy.tab' 

#If the assembly fasta is not already in the output dir, we need to copy it there so that the docker container can see it
output_dir_absolute = os.path.abspath(output_dir)
fasta_assembly_absolute = os.path.abspath(fasta_assembly)
fasta_directory = '/'.join(fasta_assembly_absolute.split('/')[:-1])
fasta_filename = fasta_assembly_absolute.split('/')[-1]

if output_dir_absolute != fasta_directory:
	# This means we need to copy the fasta assembly to the output directory
	run_command('cp ' + fasta_assembly_absolute + ' ' + output_dir_absolute + '/')

#If a taxonomy table is given, and it already exists, we will need to copy it to the output directory if it isn't already there
if taxonomy_table_path:
	taxonomy_table_path_absolute = os.path.abspath(taxonomy_table_path)
	taxonomy_directory = '/'.join(taxonomy_table_path_absolute.split('/')[:-1])
	taxonomy_table_filename = taxonomy_table_path_absolute.split('/')[-1]
	if os.path.isfile(taxonomy_table_path):
		if output_dir_absolute != taxonomy_directory:
			run_command('cp ' + taxonomy_table_path_absolute + ' ' + output_dir_absolute + '/')

#If a coverage table is give, it must already exist.
#If it exists, then it also should be copied to the output directory if it isn't already there
if cov_table:
	if not os.path.isfile(cov_table):
		print('Error! Could not find coverage table at the following path: ' + cov_table)
		exit(1)
	else:
		cov_table_absolute = os.path.abspath(cov_table)
		cov_table_directory = '/'.join(cov_table_absolute.split('/')[:-1])
		cov_table_filename = cov_table_absolute.split('/')[-1]
		if output_dir_absolute != cov_table_directory:
			run_command('cp ' + cov_table_absolute + ' ' + output_dir_absolute + '/')

# Construct run_autometa.py command to pass to the docker container
autometa_command = 'run_autometa.py --assembly /output/{} --processors {} --length_cutoff {} --completeness_cutoff {} --kingdom {} --output_dir /output'.format(\
	fasta_filename, processors, length_cutoff, cluster_completeness, kingdom)

# Add optional arguments
if taxonomy_table_path:
	autometa_command = autometa_command + ' --taxonomy_table /output/{}'.format(taxonomy_table_filename)

if db_dir_path:
	autometa_command = autometa_command + ' --db_dir /databases'

if do_ML_recruitment:
	autometa_command = autometa_command + ' --ML_recruitment'

if make_tax_table:
	autometa_command = autometa_command + ' --maketaxtable'

if cov_table:
	#Note: if this is specified, we've already checked it exists and is now in the output dir
	autometa_command = autometa_command + ' --cov_table /output/{}'.format(cov_table_filename)

# Construct Docker run command
docker_command = 'docker run --volume {}:/output:rw --detach=false --rm'.format(output_dir_absolute)

if db_dir_path:
	docker_command = docker_command + ' --volume {}:/databases:rw'.format(db_dir_path_absolute)

docker_command = docker_command + ' autometa:run'

# Incorporate autometa command
docker_command = docker_command + ' ' + autometa_command

# Run docker container
run_command(docker_command)
