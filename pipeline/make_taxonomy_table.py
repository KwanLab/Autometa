#!/usr/bin/env python

import subprocess
import argparse
import os
import pandas as pd
from Bio import SeqIO
import sys

def readable_dir(prospective_dir):
  if not os.path.isdir(prospective_dir):
    raise Exception("readable_dir: \"{0}\" is not a valid path".format(prospective_dir))
  if os.access(prospective_dir, os.R_OK):
    return prospective_dir
  else:
    raise Exception("readable_dir: \"{0}\" is not a readable directory".format(prospective_dir))

#argument parser
parser = argparse.ArgumentParser(description="Script to generate the contig taxonomy table.", epilog="Output will be directed to recursive_dbscan.py")
parser.add_argument('-a', metavar='assembly', help='assembly.fasta', required=True)
parser.add_argument('-p', metavar='processors', help='Processors to use. If not specified, will infer', default=1)
parser.add_argument('-n', metavar='NR Diamond db', help='Diamond formatted non-redundant (NR) protein database', default="/media/box2/nr_old/nr") #Need to update default
parser.add_argument('-t', metavar='database directory', help='Path to directory with taxdump files.', required=True, type=readable_dir)
parser.add_argument('-l', metavar='cutoff length', help='Contig length cutoff to consider for binning.\
 Default is 10,000 bp.', default=10000, type = int)
parser.add_argument("-update", required=False, action='store_true', help='Adds/Updates nodes.dmp, names.dmp and accession2taxid files within taxdump directory specified with \"-t\"')

args = vars(parser.parse_args())

num_processors = args['p']
length_cutoff = args['l']
fasta_path = args['a']
fasta_assembly_prefix = os.path.splitext(os.path.basename(args['a']))[0]

def length_trim(fasta_path,fasta_prefix,length_cutoff):
	#Trim the length of fasta file
	outfile_name = str(fasta_prefix) + "_filtered.fasta"
	subprocess.call("{}/fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta_path, length_cutoff, outfile_name), shell = True)
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
    subprocess.call("diamond blastp --query {}.faa --db {} --evalue 1e-5 --max-target-seqs 200 -p {} --daa {}".format(prodigal_output, diamond_database_path, num_processors, prodigal_daa), shell = True)
    subprocess.call("diamond view -a {} -f tab -o {}".format(prodigal_daa, view_output), shell = True)
    return view_output

	#return  view_output
	#might want to change name of outputfile
"""
#blast2lca using accession numbers#
def run_blast2lca(input_file, taxdump_path):
	output = input_file.rstrip(".tab") + ".lca"
	if os.path.isfile(output):
		print "{} file already exists!".format(output)
		print "Continuing to next step..."
	else:
		subprocess.call("{}/lca.py database_directory {} {} > {}".format(pipeline_path, taxdump_dir_path, input_file, output), shell = True)
	return output

"""

def run_blast2lca(input_file,taxdump_path):
	output = input_file.rstrip(".tab") + ".lca"
	if os.path.isfile(output):
		print "{} file already exists!".format(output)
		print "Continuing to next step..."
	else:
		subprocess.call("blast2lca -savemem -dict {}/gi_taxid.bin -nodes {}/nodes.dmp -names {}/names.dmp {} > {}"\
			.format(taxdump_path,input_file, output), shell = True)
	return output

def run_taxonomy(pipeline_path, assembly_path, tax_table_path, taxdump_dir_path): #Have to update this
	initial_table_path = assembly_path + '.tab'
	subprocess.call("{}/make_contig_table.py {} {}".format(pipeline_path, assembly_path, initial_table_path), shell = True)
	subprocess.call("{}/add_contig_taxonomy.py {} {} {} taxonomy.tab".format(pipeline_path, initial_table_path, tax_table_path, taxdump_dir_path), shell = True)
	#contig_table_path, tax_table_path, taxdump_dir_path, output_file_path
	return 'taxonomy.tab'

#diamond_path = subprocess.check_output('find ~ -name "diamond"', shell=True).rstrip("\n") # assume diamond is in the path
taxdump_dir_path = args['t']#'/home/jkwan/blast2lca_taxdb'
prodigal_output = fasta_assembly_prefix + "_filtered.orfs"
prodigal_daa = prodigal_output + ".daa"
pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)
#diamond_database_path = subprocess.check_output('find /mnt/not_backed_up/nr_diamond/ -name "nr.dmnd"', shell=True).strip("\n")
diamond_database_path = args['n'] #/media/box2/nr_old/nr'
#add_contig_path = pipeline_path
filtered_assembly = fasta_assembly_prefix + "_filtered.fasta"
names_dmp_path = taxdump_dir_path + '/names.dmp'
nodes_dmp_path = taxdump_dir_path + '/nodes.dmp'
accession2taxid_path = taxdump_dir_path + '/prot.accession2taxid'
taxdump_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
accession2taxid_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"

if args['update']:
    os.system('wget -v %s -O %s/taxdump.tar.gz' % (taxdump_url, taxdump_dir_path))
    os.system('tar -xvzf %s/taxdump.tar.gz -C %s names.dmp nodes.dmp' % (taxdump_dir_path, taxdump_dir_path))
    os.system('rm %s/taxdum.tar.gz' % taxdump_dir_path)
    os.system('wget -v %s -O %s.gz' % (accession2taxid_url, accession2taxid_path))
    print("Gunzipping prot.accession2taxid gzipped file\nThis may take some time...")
    os.system('gunzip -9vNf %s.gz > %s' % (accession2taxid_path, accession2taxid_path))
    print("\nnodes.dmp, names.dmp and prot.accession2taxid files updated in %s\nResuming..." % taxdump_dir_path)

if not os.path.isfile(prodigal_output + ".faa"):
	print "Prodigal output not found. Running prodigal..."
	#Check for file and if it doesn't exist run make_marker_table
	length_trim(fasta_path, fasta_assembly_prefix, length_cutoff)
	run_prodigal(filtered_assembly)

if not os.path.isfile(prodigal_output + ".daa"):
	print "Could not find {}. Running diamond blast... ".format(prodigal_output + ".daa")
	#diamond_output =
	diamond_output = run_diamond(prodigal_output, diamond_database_path, num_processors, prodigal_daa)
elif os.stat(prodigal_output + ".daa").st_size == 0:
	print "{} files empty. Re-running diamond blast...".format(prodigal_output + ".tab")
	diamond_output = run_diamond(prodigal_output, diamond_database_path, num_processors, prodigal_daa)
else:
	diamond_output = prodigal_output + ".tab"

if not os.path.isfile(prodigal_output + ".lca"):
    blast2lca_output = run_blast2lca(diamond_output,taxdump_dir_path)
else:
    blast2lca_output = prodigal_output + ".lca"

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
