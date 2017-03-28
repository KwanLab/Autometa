#!/usr/bin/env python

import subprocess
import argparse
import os
import pandas as pd
from Bio import SeqIO
import sys

#argument parser
parser = argparse.ArgumentParser(description="Script to generate the contig taxonomy table.")
parser.add_argument('-a','--assembly', help='assembly.fasta', required=True)
parser.add_argument('-p','--processors', help='assembly.fasta', default=1)
parser.add_argument('-l','--length_cutoff', help='Contig length cutoff to consider for binning.\
 Default is 10,000 bp.', default=10000, type = int)
args = vars(parser.parse_args())

num_processors = args['processors']
length_cutoff = args['length_cutoff']
fasta_assembly = os.path.basename(args['assembly'])


def length_trim(fasta,length_cutoff):
	#Trim the length of fasta file
	outfile_name = str(args['assembly'].split("/")[-1].split(".")[0]) + "_filtered.fasta"
	subprocess.call("{}/fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff, outfile_name), shell = True)
	return outfile_name

def run_prodigal(path_to_assembly):
	#When "shell = True", need to give one string, not a list
	prodigal_output = path_to_assembly.split('.')[0] + '.orfs.faa'
	if os.path.isfile(prodigal_output):
		print "{} file already exists!".format(prodigal_output)
		print "Continuing to next step..."
	else:
		subprocess.call(" ".join(['prodigal ','-i ' + path_to_assembly, '-a ' + path_to_assembly.split(".")[0] +\
	 	'.orfs.faa','-p meta', '-m', '-o ' + path_to_assembly.split(".")[0] + '.txt']), shell = True)

def run_diamond(prodigal_output, diamond_database_path, num_processors, prodigal_daa):
	view_output = prodigal_output + ".tab"
	current_dir = os.getcwd()
	tmp_dir_path = current_dir + '/tmp'
	if not os.path.isdir(tmp_dir_path):
		os.makedirs(tmp_dir_path) # This will give an error if the path exists but is a file instead of a dir
	if os.path.isfile(prodigal_daa):
		print "{} file already exists!".format(prodigal_daa)
		print 'Continuing to next step...'
	else:
		subprocess.call("diamond blastp --query {}.faa --db {} --evalue 1e-5 --max-target-seqs 200 -p {} --daa {} -t {}"\
			.format(prodigal_output, diamond_database_path, num_processors, prodigal_daa, tmp_dir_path), shell = True)

	if os.path.isfile(view_output):
		print "{} file already exists!".format(view_output)
		print "Continuing to next setp..."
	else:
		subprocess.call("diamond view -a {} -f tab -o {}".format(prodigal_daa, view_output), shell = True)
	return view_output
	#return  view_output
	#might want to change name of outputfile

def run_blast2lca(input_file):
	output = input_file.rstrip(".tab") + ".lca"
	if os.path.isfile(output):
		print "{} file already exists!".format(output)
		print "Continuing to next step..."
	else:
		subprocess.call("/home/ijmiller/go/bin/blast2lca -savemem -dict /home/jkwan/blast2lca_taxdb/gi_taxid.bin -nodes /home/jkwan/blast2lca_taxdb/nodes.dmp -names /home/jkwan/blast2lca_taxdb/names.dmp {} > {}"\
			.format(input_file, output), shell = True)
	return output

def run_taxonomy(pipeline_path, assembly_path, tax_table_path, taxdump_dir_path): #Have to update this
	initial_table_path = assembly_path + '.tab'
	subprocess.call("{}/make_contig_table.py {} {}".format(pipeline_path, assembly_path, initial_table_path), shell = True)
	subprocess.call("{}/add_contig_taxonomy.py {} {} {} taxonomy.tab".format(pipeline_path, initial_table_path ,tax_table_path, taxdump_dir_path), shell = True)
	#contig_table_path, tax_table_path, taxdump_dir_path, output_file_path
	return 'taxonomy.tab'

#diamond_path = subprocess.check_output('find ~ -name "diamond"', shell=True).rstrip("\n") # assume diamond is in the path
taxdump_dir_path = '/home/jkwan/blast2lca_taxdb'
prodigal_output = args['assembly'].split('.')[0] + "_filtered.orfs"
prodigal_daa = prodigal_output + ".daa"
pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)
#diamond_database_path = subprocess.check_output('find /mnt/not_backed_up/nr_diamond/ -name "nr.dmnd"', shell=True).strip("\n")
diamond_database_path = '/mnt/not_backed_up/nr_diamond/nr'
#add_contig_path = pipeline_path
filtered_assembly = str(args['assembly'].split("/")[-1].split(".")[0]) + "_filtered.fasta"

if not os.path.isfile(prodigal_output + ".faa"):
	print "Prodigal output not found. Running prodigal..."
	#Check for file and if it doesn't exist run make_marker_table
	length_trim(fasta_assembly, length_cutoff)
	run_prodigal(filtered_assembly)


if not os.path.isfile(prodigal_output + ".tab"):
	print "Running diamond blast... "
	#diamond_output =
	diamond_output = run_diamond(prodigal_output, diamond_database_path, num_processors, prodigal_daa)
else:
	diamond_output = prodigal_output + ".tab"

blast2lca_output = run_blast2lca(diamond_output)
print "Running add_contig_taxonomy.py... "
taxonomy_table = run_taxonomy(pipeline_path, filtered_assembly, blast2lca_output, taxdump_dir_path)

# Split the original contigs into sets for each kingdom
taxonomy_pd = pd.read_table(taxonomy_table)
categorized_seq_objects = {}
all_seq_records = {}

# Load fasta file
for seq_record in SeqIO.parse(filtered_assembly, 'fasta'):
	all_seq_records[seq_record.id] = seq_record

for i, row in taxonomy_pd.iterrows():
	kingdom = row['kingdom']
	contig = row['contig']
	if kingdom in categorized_seq_objects:
		categorized_seq_objects[kingdom].append(all_seq_records[contig])
	else:
		categorized_seq_objects[kingdom] = [ all_seq_records[contig] ]

# Now we write the component fasta files
for kingdom in categorized_seq_objects:
	seq_list = categorized_seq_objects[kingdom]
	output_path = kingdom + '.fasta'
	SeqIO.write(seq_list, output_path, 'fasta')

print "Done!"
