#!/usr/bin/env python

# Program to proprocess contigs into kingdoms - Bacteria, Archeaea, Eukaryota and unclassified.
# This is an optional step to be carried out before run_autometa.py.  Recommended for host-associated metagenomes.

import sys
import os
import pandas as pd
import csv

input_fasta = sys.argv[1]

def length_trim(fasta,length_cutoff):
	#will need to update path of this perl script
	outfile_name = str(args['assembly'].split("/")[-1].split(".")[0]) + "_filtered.fasta"
	subprocess.call("{}/fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff,outfile_name), shell = True)
	return outfile_name

def make_contig_table(fasta):
	#looks like this script is assuming contigs from a spades assembly
	output_table_name = str(fasta).split('.')[0] + ".tab"
	subprocess.call("{}/make_contig_table.py {} {}".format(pipeline_path,fasta,output_table_name), shell = True)
	return output_table_name

def run_prodigal(path_to_assembly):
	#When "shell = True", need to give one string, not a list
	subprocess.call(" ".join(['prodigal','-i ' + path_to_assembly, '-a ' + path_to_assembly.split(".")[0] +\
	 '.orfs.faa','-p meta', '-m', '-o ' + path_to_assembly.split(".")[0] + '.txt']), shell = True)

pipeline_path = sys.path[0]

filtered_assembly = length_trim(input_fasta,3000)
contig_table = make_contig_table(filtered_assembly)

# Run prodigal
run_prodigal(filtered_assembly)
