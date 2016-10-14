#!/usr/bin/env python

import argparse
import subprocess
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
parser = argparse.ArgumentParser(description='Script to tabulate paired-end read \
coverage using samtools and Bedtools. Note: it is recommended that you quality \
filter your *.fastq files prior to calculating read coverage.')
parser.add_argument('-a','--assembly', help='Input assembly file', required=True)
parser.add_argument('-p','--processors', help='Number of processors', default=1)
parser.add_argument('-F','--forward_reads', help='Paired (*.fastq) forward reads', required=True)
parser.add_argument('-R','--reverse_reads', help='Paired (*.fastq) reverse reads', required=True)
parser.add_argument('-o','--out', help='Tab delimited table for contig and average\
    read coverage', default="outfile.tab")
args = vars(parser.parse_args())

def run_bowtie2(path_to_assembly,path_to_F_reads,path_to_R_reads,num_processors=1):
	#When "shell = True", need to give one string, not a list
    #Build bowtie2 database
    bowtie2_db = path_to_assembly.split(".")[0]
    subprocess.call(" ".join(['bowtie2-build ', path_to_assembly, bowtie2_db]), shell = True)
    #Run bowtie2 alignment
    sam_file_name = bowtie2_db + '.sam'
    subprocess.call(" ".join(['bowtie2','-x ' + bowtie2_db, '-1 ' + path_to_F_reads,\
       '-2 ' + path_to_R_reads, '-q','--phred33','--very-sensitive',\
       '--no-unal','-p ',str(num_processors),'-S ' + sam_file_name]), shell = True)
    return sam_file_name

#Check for dependencies in $PATH

#1. Align the comparison dataset reads to your assembly with bowtie2.
assembly_file = args['assembly']

sam_file = run_bowtie2(args['assembly'],args['forward_reads'],\
    args['reverse_reads'],args['processors'])

#2. Convert the SAM file to a sorted BAM file, and create a BAM index.
sorted_bam_file = assembly_file.split(".")[0] + ".sort"
subprocess.call(" ".join(['samtools view ','-bS ' + sam_file, ' | ', \
    'samtools sort ', '-o ', sorted_bam_file]), shell = True)

#3. Tabulate the average coverage of each contig.

#Run fasta_length_table.pl
contig_length_tab_file = assembly_file.split(".")[0] + ".tab"
subprocess.call(" ".join(['fasta_length_table.pl ', assembly_file, '> ',\
    contig_length_tab_file]), shell = True)

#Run genomeCoverageBed
genome_bed_file = assembly_file.split(".")[0] + ".txt"
subprocess.call(" ".join(['genomeCoverageBed ', '-ibam ', sorted_bam_file, \
    '-g ' + contig_length_tab_file, '> ' + genome_bed_file]), shell = True)

#Build final table
outfile = args['out']
subprocess.call(" ".join(['contig_coverage_from_bedtools.pl ', genome_bed_file, \
    '> ' + outfile]), shell = True)

#Future features: check for .fastq.gz, dependencies in path, default processors
