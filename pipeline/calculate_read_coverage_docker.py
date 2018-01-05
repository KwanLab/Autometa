#!/usr/bin/env python

# Docker wrapper for calculate_read_coverage.py
# Aims to be invisible to the user (use in exactly the same way as calculate_read_coverage.py)

import argparse
import os
import subprocess

def run_command(command_string, stdout_path = None):
	# Function that checks if a command ran properly. If it didn't, then print an error message then quit
	if stdout_path:
		f = open(stdout_path, 'w')
		exit_code = subprocess.call(command_string, stdout=f, shell=True)
		f.close()
	else:
		exit_code = subprocess.call(command_string, shell=True)

	if exit_code != 0:
		print('calculate_read_coverage_docker.py: Error, the command:')
		print(command_string)
		print('failed, with exit code ' + str(exit_code))
		exit(1)

#argument parser
parser = argparse.ArgumentParser(description='Script to tabulate paired-end or single read \
coverage using Samtools and Bedtools. Note: it is recommended that you quality \
filter your *.fastq or *fastq.gz files prior to calculating read coverage.')
parser.add_argument('-a','--assembly', metavar='<assembly.fasta>', help='Path to input assembly file', required=True)
parser.add_argument('-p','--processors', metavar='<int>', help='Number of processors', type=int, default=1)
parser.add_argument('-F','--forward_reads', metavar='<reads.fastq|reads.fastq.gz>', help='Paired (*.fastq|*.fastq.gz) forward reads (must be in same order as reverse list)', nargs='*')
parser.add_argument('-R','--reverse_reads', metavar='<reads.fastq|reads.fastq.gz>', help='Paired (*.fastq|*.fastq.gz) reverse reads (must be in same order as forward list)', nargs='*')
parser.add_argument('-S','--single_reads', metavar='<reads.fastq|reads.fastq.gz>', help='Single (*.fastq|*.fastq.gz) reads', nargs='*')
parser.add_argument('-o','--out', metavar='<output.tab>', help='Tab delimited table for contig and average\
    read coverage', default="outfile.tab")
args = vars(parser.parse_args())

assembly_file = args['assembly']
processors = args['processors']
forward_read_path_list = args['forward_reads']
reverse_read_path_list = args['reverse_reads']
single_read_path_list = args['single_reads']
outfile = args['out']

if not os.path.isfile(assembly_file):
	print ('Error! Could not find assembly file at the following path: ' + assembly_file)
	exit(1)

# Check that at least one of the above lists is not empty
if not (forward_read_path_list or reverse_read_path_list or single_read_path_list):
    print('Error! You need to specify some reads, with -F/-R (paired) and/or -S (single)')
    exit(1)

# Check that all read paths exist
concatenated_read_path_list = forward_read_path_list + reverse_read_path_list + single_read_path_list
for path in concatenated_read_path_list:
    if not os.path.isfile(path):
        print('Error! Cannot find the read file: ' + path)
        exit(1)

# We will use the path for the output.tab as showing us which output dir to use
outfile_absolute = os.path.abspath(outfile)
output_dir = '/'.join(outfile_absolute.split('/')[:-1])
outfile_filename = outfile_absolute.split('/')[-1]

# We have to create the output dir if it doesn't exist
if not os.path.isdir(output_dir):
	os.makedirs(output_dir)

# If the fasta file is not in the output dir, we copy it to the output dir
assembly_file_absolute = os.path.abspath(assembly_file)
assembly_dir = '/'.join(assembly_file_absolute.split('/')[:-1])
assembly_filename = assembly_file_absolute.split('/')[-1]

if output_dir != assembly_dir:
	run_command('cp ' + assembly_file_absolute + ' ' + output_dir + '/')

# For the reads, we need to make sure they all exist and that we hook up the docker container to the correct directories
read_directories = dict()
read_directory_counter = 0
docker_forward_read_path_list = list()
docker_reverse_read_path_list = list()
docker_single_read_path_list = list()

if forward_read_path_list:
	for path in forward_read_path_list:
		read_abs_path = os.path.abspath(path)
		read_dir = '/'.join(read_abs_path.split('/')[:-1])
		read_filename = read_abs_path.split('/')[-1]
		if read_dir not in read_directories:
			read_directory_counter += 1
			docker_dir_name = 'reads_dir_' + str(read_directory_counter)
			read_directories[read_dir] = docker_dir_name
		docker_forward_read_path_list.append(read_directories[read_dir] + '/' + read_filename)
	for path in reverse_read_path_list:
		read_abs_path = os.path.abspath(path)
		read_dir = '/'.join(read_abs_path.split('/')[:-1])
		read_filename = read_abs_path.split('/')[-1]
		if read_dir not in read_directories:
			read_directory_counter += 1
			docker_dir_name = 'reads_dir_' + str(read_directory_counter)
			read_directories[read_dir] = docker_dir_name
		docker_reverse_read_path_list.append(read_directories[read_dir] + '/' + read_filename)

if single_read_path_list:
	for path in single_read_path_list:
		read_abs_path = os.path.abspath(path)
		read_dir = '/'.join(read_abs_path.split('/')[:-1])
		read_filename = read_abs_path.split('/')[-1]
		if read_dir not in read_directories:
			read_directory_counter += 1
			docker_dir_name = 'reads_dir_' + str(read_directory_counter)
			read_directories[read_dir] = docker_dir_name
		docker_single_read_path_list.append(read_directories[read_dir] + '/' + read_filename)

# Construct the calculate_read_coverage.py command
calculate_read_coverage_command = 'calculate_read_coverage.py --assembly /output/{} --processors {} --out /output/{}'.format(\
	assembly_filename, processors, outfile_filename)

if docker_forward_read_path_list:
	calculate_read_coverage_command = calculate_read_coverage_command + ' --forward_reads ' + ' '.join(docker_forward_read_path_list)
	calculate_read_coverage_command = calculate_read_coverage_command + ' --reverse_reads ' + ' '.join(docker_reverse_read_path_list)

if docker_single_read_path_list:
	calculate_read_coverage_command = calculate_read_coverage_command + ' --single_reads ' + ' '.join(docker_single_read_path_list)

# Construct the docker run command
docker_command = 'docker run --volume {}:/output:rw --detach=false --rm'.format(output_dir)

# Now we need to add how every many read dirs there are
for host_path in read_directories:
	docker_command = docker_command + ' --volume {}:{}:rw'.format(host_path, read_directories[host_path])

docker_command = docker_command + ' autometa:run'

# Incorporate autometa command
docker_command = docker_command + ' ' + calculate_read_coverage_command

# Run docker container
run_command(docker_command)