#!/usr/bin/env python

import argparse
import subprocess
import os
#import multiprocessing

#Input: fastq (or fastq.gz?) reads, metagenomic assembly contigs as fasta file
#Output: "contig\tread_coverage"

#See: /home/ijmiller/slurm_scripts/2016/October/7-Oct-16_Bpac5Tips_coverage_calculation.sh
#Dependencies: bowtie2 (2.2.5), samtools (0.1.18), \
    #Bedtools (genomeCoverageBed: v2.17.0), contig_coverage_from_bedtools.pl,\
    #fasta_length_table.pl, python modules: suprocess, argparse

#1. Align the reads to your assembly with bowtie2.
#2. Convert the SAM file to a sorted BAM file, and create a BAM index.
#3. Tabulate the average coverage of each contig.

#num_cpus = multiprocessing.cpu_count()
#use_cpus = int(round((num_cpus / 2)))

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

def run_command(command_string, stdout_path = None):
    # Function that checks if a command ran properly. If it didn't, then print an error message then quit
    if stdout_path:
        f = open(stdout_path, 'w')
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print('calculate_read_coverage.py: Error, the command:')
        print(command_string)
        print('failed, with exit code ' + str(exit_code))
        exit(1)

def run_bowtie2(path_to_assembly,paths_to_F_reads,paths_to_R_reads,paths_to_S_reads,num_processors=1):
    # Check the number/type of reads we have as input
    if not (len(paths_to_F_reads) or len(paths_to_R_reads) or len(paths_to_S_reads)):
        print ('Error, you need to specify at least one read file')
        quit()

    if not (len(paths_to_F_reads) == len(paths_to_R_reads)):
        print('Error, you need to specify the same number of forward and reverse read files')
        quit()

    if paths_to_F_reads:
        F_read_path_string = ','.join(paths_to_F_reads)
        R_read_path_string = ','.join(paths_to_R_reads)

    if paths_to_S_reads:
        S_read_path_string = ','.join(paths_to_S_reads)

	#When "shell = True", need to give one string, not a list
    #Build bowtie2 database
    bowtie2_db = os.path.splitext(os.path.basename(path_to_assembly))[0]
    run_command(" ".join(['bowtie2-build ',path_to_assembly, bowtie2_db]))
    #Run bowtie2 alignment
    sam_file_name = bowtie2_db + '.sam'

    bowtie2_command_string = 'bowtie2 -x ' + bowtie2_db + ' -q --phred33 --very-sensitive --no-unal -p ' + str(num_processors) + ' -S ' + sam_file_name

    if paths_to_F_reads:
        bowtie2_command_string = bowtie2_command_string + ' -1 ' + F_read_path_string
        bowtie2_command_string = bowtie2_command_string + ' -2 ' + R_read_path_string

    if paths_to_S_reads:
        bowtie2_command_string = bowtie2_command_string + ' -U ' + S_read_path_string

    run_command(bowtie2_command_string)
    return sam_file_name

#Check for dependencies in $PATH

#1. Align the comparison dataset reads to your assembly with bowtie2.
assembly_file = args['assembly']
assembly_file_prefix = os.path.splitext(os.path.basename(assembly_file))[0]

forward_read_path_list = args['forward_reads']
reverse_read_path_list = args['reverse_reads']
single_read_path_list = args['single_reads']

if not os.path.isfile(assembly_file):
    print ('Error! Could not find assembly file at the following path: ' + assembly_file)
    exit(1)

# Check that at least one of the above lists is not empty
if not (forward_read_path_list or reverse_read_path_list or single_read_path_list):
    print('Error! You need to specify some reads, with -F/-R (paired) and/or -S (single)')
    exit(1)

# Check that all read paths exist
if forward_read_path_list:
    concatenated_read_path_list = forward_read_path_list
if reverse_read_path_list:
    concatenated_read_path_list = concatenated_read_path_list + reverse_read_path_list
if single_read_path_list:
    concatenated_read_path_list = concatenated_read_path_list + single_read_path_list
for path in concatenated_read_path_list:
    if not os.path.isfile(path):
        print('Error! Cannot find the read file: ' + path)
        exit(1)

sam_file = run_bowtie2(assembly_file,forward_read_path_list,reverse_read_path_list,single_read_path_list,args['processors'])

#2. Convert the SAM file to a sorted BAM file, and create a BAM index.
sorted_bam_file = assembly_file_prefix + ".sort.bam"
run_command('samtools view -bS {} | samtools sort -o {}'.format(sam_file, sorted_bam_file))

#Clean up the SAM file, which will be a lot larger than the sorted BAM file
run_command('rm ' + sam_file)

#3. Tabulate the average coverage of each contig.

#Run fasta_length_table.pl
contig_length_tab_file = assembly_file_prefix + ".tab"
run_command('fasta_length_table.pl ' + assembly_file, contig_length_tab_file)

#Run genomeCoverageBed
genome_bed_file = assembly_file_prefix + ".txt"
run_command('genomeCoverageBed -ibam {} -g {}'.format(sorted_bam_file, contig_length_tab_file), genome_bed_file)

#Build final table
outfile = args['out']
run_command('contig_coverage_from_bedtools.pl ' + genome_bed_file, outfile)

#Future features: check for .fastq.gz, dependencies in path, default processors
