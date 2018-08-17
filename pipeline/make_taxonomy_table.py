#!/usr/bin/env python

# Copyright 2018 Ian J. Miller, Evan R. Rees, Izaak Miller, Jason C. Kwan
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

import subprocess
import argparse
import os
import pandas as pd
from Bio import SeqIO
import sys
import urllib2

def run_command(command_string, stdout_path = None):
	# Function that checks if a command ran properly. If it didn't, then print an error message then quit
	print('make_taxonomy_table.py, run_command: ' + command_string)
	if stdout_path:
		f = open(stdout_path, 'w')
		exit_code = subprocess.call(command_string, stdout=f, shell=True)
		f.close()
	else:
		exit_code = subprocess.call(command_string, shell=True)

	if exit_code != 0:
		print('Make_taxonomy_table.py: Error, the command:')
		print(command_string)
		print('failed, with exit code ' + str(exit_code))
		exit(1)

def run_command_return(command_string, stdout_path = None):
	# Function that checks if a command ran properly. If it didn't, then print an error message then quit
	print('make_taxonomy_table.py, run_command: ' + command_string)
	if stdout_path:
		f = open(stdout_path, 'w')
		exit_code = subprocess.call(command_string, stdout=f, shell=True)
		f.close()
	else:
		exit_code = subprocess.call(command_string, shell=True)

	return exit_code

def cythonize_lca_functions():
	print("{}/lca_functions.so not found, cythonizing lca_function.pyx for make_taxonomy_table.py".format(pipeline_path))
	current_dir = os.getcwd()
	os.chdir(pipeline_path)
	run_command("python setup_lca_functions.py build_ext --inplace")
	os.chdir(current_dir)

def download_file(destination_dir, file_url, md5_url):
	filename = file_url.split('/')[-1]
	md5name = md5_url.split('/')[-1]

	md5check = False

	while md5check == False:
		run_command('wget {} -O {}'.format(file_url, destination_dir + '/' + filename))
		run_command('wget {} -O {}'.format(md5_url, destination_dir + '/' + md5name))

		downloaded_md5 = subprocess.check_output(['md5sum', destination_dir + '/' + filename]).split(' ')[0]

		with open(destination_dir + '/' + md5name, 'r') as check_md5_file:
			check_md5 = check_md5_file.readline().split(' ')[0]

		if downloaded_md5 == check_md5:
			md5check = True

def md5IsCurrent(local_md5_path, remote_md5_url):
	remote_md5_handle = urllib2.urlopen(remote_md5_url)
	remote_md5 = remote_md5_handle.readline().split(' ')[0]

	with open(local_md5_path, 'r') as local_md5_file:
		local_md5 = local_md5_file.readline().split(' ')[0]

	if local_md5 == remote_md5:
		return True
	else:
		return False

def update_dbs(database_path, db='all'):
	"""Updates databases for AutoMeta usage"""
	taxdump_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
	taxdump_md5_url = taxdump_url+".md5"
	accession2taxid_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
	accession2taxid_md5_url = accession2taxid_url+".md5"
	nr_db_url = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
	nr_db_md5_url = nr_db_url+".md5"
	#Downloading files for db population
	if db == 'all' or db == 'nr':
		# First download nr if we don't yet have it OR it is not up to date
		if os.path.isfile(database_path + '/nr.gz.md5'):
			if not md5IsCurrent(database_path + '/nr.gz.md5', nr_db_md5_url):
				print("updating nr.dmnd")
				download_file(database_path, nr_db_url, nr_db_md5_url)
		else:
			print("updating nr.dmnd")
			download_file(database_path, nr_db_url, nr_db_md5_url)

		# Now we make the diamond database
		print("building nr.dmnd database, this may take some time")
		returnCode = subprocess.call("diamond makedb --in {} --db {}/nr -p {}".format(database_path+'/nr.gz', database_path, num_processors), shell = True)
		if returnCode == 0: # i.e. job was successful
		#Make an md5 file to signal that we have built the database successfully
			run_command('rm {}/nr.gz'.format(database_path))
			print("nr.dmnd updated")
	if db == 'all' or db == 'acc2taxid':
		# Download prot.accession2taxid.gz only if the version we have is not current
		if os.path.isfile(database_path + '/prot.accession2taxid.gz.md5'):
			if not md5IsCurrent(database_path + '/prot.accession2taxid.gz.md5', accession2taxid_md5_url):
				print("updating prot.accession2taxid")
				download_file(database_path, accession2taxid_url, accession2taxid_md5_url)
		else:
			print("updating prot.accession2taxid")
			download_file(database_path, accession2taxid_url, accession2taxid_md5_url)

		if os.path.isfile(database_path + '/prot.accession2taxid.gz'):
			print("Gunzipping prot.accession2taxid gzipped file\nThis may take some time...")
			run_command('gunzip -9vNf {}'.format(database_path + '/prot.accession2taxid.gz'), database_path + '/prot.accession2taxid')
			print("prot.accession2taxid updated")
	if db == 'all' or db == 'taxdump':
		# Download taxdump only if the version we have is not current
		if os.path.isfile(database_path + '/taxdump.tar.gz.md5'):
			if not md5IsCurrent(database_path + '/taxdump.tar.gz.md5', taxdump_md5_url):
				print("updating nodes.dmp and names.dmp")
				download_file(database_path, taxdump_url, taxdump_md5_url)
		else:
			print("updating nodes.dmp and names.dmp")
			download_file(database_path, taxdump_url, taxdump_md5_url)

		if os.path.isfile(database_path + '/taxdump.tar.gz'):
			run_command('tar -xzf {}/taxdump.tar.gz -C {} names.dmp nodes.dmp'.format(database_path, database_path))
			run_command('rm {}/taxdump.tar.gz'.format(database_path))
			print("nodes.dmp and names.dmp updated")

def length_trim(fasta_path,length_cutoff):
	input_filename = '.'.join(os.path.abspath(fasta_path).split('/')[-1].split('.')[:-1])
	#Trim the length of fasta file
	outfile_path = output_dir + '/' + input_filename + "_filtered.fasta"
	run_command("{}/fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta_path, length_cutoff, outfile_path))
	return outfile_path

def run_prodigal(path_to_assembly):
	assembly_filename = os.path.abspath(path_to_assembly).split('/')[-1]
	assembly_basename = '.'.join(assembly_filename.split('.')[:-1])
	output_path = output_dir + '/' + assembly_basename + '.orfs.faa'
	if os.path.isfile(output_path):
		print "{} file already exists!".format(output_path)
		print "Continuing to next step..."
	else:
		run_command('prodigal -i {} -a {}/{}.orfs.faa -p meta -m -o {}/{}.txt'.format(path_to_assembly, output_dir, assembly_basename, output_dir, assembly_basename))

def run_diamond(prodigal_output, diamond_db_path, num_processors, prodigal_daa):
	view_output = prodigal_output + ".tab"
	current_dir = os.getcwd()
	tmp_dir_path = current_dir + '/tmp'
	if not os.path.isdir(tmp_dir_path):
		os.makedirs(tmp_dir_path) # This will give an error if the path exists but is a file instead of a dir
	error = run_command_return("diamond blastp --query {}.faa --db {} --evalue 1e-5 --max-target-seqs 200 -p {} --daa {} -t {}".format(prodigal_output, diamond_db_path, num_processors, prodigal_daa,tmp_dir_path))
	# If there is an error, attempt to rebuild NR
	if error:
		update_dbs(db_dir_path, 'nr')
		run_command("diamond blastp --query {}.faa --db {} --evalue 1e-5 --max-target-seqs 200 -p {} --daa {} -t {}".format(prodigal_output, diamond_db_path, num_processors, prodigal_daa,tmp_dir_path))

	run_command("diamond view -a {} -f tab -o {}".format(prodigal_daa, view_output))
	return view_output

#blast2lca using accession numbers#
def run_blast2lca(input_file, taxdump_path):
	output = output_dir + '/' + '.'.join(os.path.abspath(input_file).split('/')[-1].split('.')[:-1]) + ".lca"
	if os.path.isfile(output) and not os.stat(prodigal_output + ".lca").st_size == 0:
		print "{} file already exists!".format(output)
		print "Continuing to next step..."
	else:
		run_command("{}/lca.py database_directory {} {}".format(pipeline_path, db_dir_path, input_file))
	return output

def run_taxonomy(pipeline_path, assembly_path, tax_table_path, db_dir_path,coverage_table): #Have to update this
	assembly_filename = assembly_path.split('/')[-1]
	initial_table_path = output_dir + '/' + assembly_filename + '.tab'

	# Only make the contig table if it doesn't already exist
	if not os.path.isfile(initial_table_path):
		if coverage_table:
			run_command("{}/make_contig_table.py -a {} -o {} -c {}".format(pipeline_path, assembly_path, initial_table_path,coverage_table))
		elif single_genome_mode:
			run_command("{}/make_contig_table.py -a {} -o {} -n".format(pipeline_path, assembly_path, initial_table_path))
		else:
			run_command("{}/make_contig_table.py -a {} -o {}".format(pipeline_path, assembly_path, initial_table_path))

	run_command("{}/add_contig_taxonomy.py {} {} {} {}/taxonomy.tab".format(pipeline_path, initial_table_path, tax_table_path, db_dir_path, output_dir))

	return output_dir + '/' + 'taxonomy.tab'

pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)

#argument parser
parser = argparse.ArgumentParser(description="Script to generate the contig taxonomy table.", epilog="Output will be directed to recursive_dbscan_output.tab")
parser.add_argument('-a', '--assembly', metavar='<assembly.fasta>', help='Path to metagenomic assembly fasta', required=True)
parser.add_argument('-p', '--processors', metavar='<int>', help='Number of processors to use.', type=int, default=1)
parser.add_argument('-db', '--db_dir', metavar='<dir>', help='Path to directory with taxdump, protein accessions and diamond (NR) protein files. If this path does not exist, will create and download files.', required=False, default=autometa_path + '/databases')
parser.add_argument('-udb', '--user_prot_db', metavar='<user_prot_db>', help='Replaces the default diamond database (nr.dmnd)', required=False)
parser.add_argument('-l', '--length_cutoff', metavar='<int>', help='Contig length cutoff to consider for binning in bp', default=10000, type = int)
parser.add_argument('-v', '--cov_table', metavar='<coverage.tab>', help="Path to coverage table made by calculate_read_coverage.py. If this is not specified then coverage information will be extracted from contig names (SPAdes format)", required=False)
parser.add_argument('-o', '--output_dir', metavar='<dir>', help='Path to directory to store output files', default='.')
parser.add_argument('-s', '--single_genome', help='Specifies single genome mode', action='store_true')
parser.add_argument('-u', '--update', required=False, action='store_true',\
 help='Checks/Adds/Updates: nodes.dmp, names.dmp, accession2taxid, nr.dmnd files within specified directory.')


args = vars(parser.parse_args())

db_dir_path = args['db_dir'].rstrip('/')
usr_prot_path = args['user_prot_db']
num_processors = args['processors']
length_cutoff = args['length_cutoff']
fasta_path = args['assembly']
cov_table = args['cov_table']
output_dir = args['output_dir']
single_genome_mode = args['single_genome']

fasta_filename = os.path.abspath(fasta_path).split('/')[-1]

prodigal_output = output_dir + '/' + '.'.join(fasta_filename.split('.')[:-1]) + "_filtered.orfs"
prodigal_daa = prodigal_output + ".daa"
filtered_assembly = output_dir + '/' + '.'.join(fasta_filename.split('.')[:-1]) + "_filtered.fasta"

# If cov_table defined, we need to check the file exists
if cov_table:
	if not os.path.isfile(cov_table):
		print("Error! Could not find coverage table at the following path: " + cov_table)
		exit(1)

# Check that output dir exists, and create it if it doesn't
if not os.path.isdir(output_dir):
	os.makedirs(output_dir)

if not os.path.isfile(pipeline_path+"/lca_functions.so"):
	cythonize_lca_functions()

if not os.path.isdir(db_dir_path):
	#Verify the 'Autometa databases' directory exists
	print('No databases directory found, creating and populating AutoMeta databases directory\nThis may take some time...')
	run_command('mkdir {}'.format(db_dir_path))
	update_dbs(db_dir_path)
elif not os.listdir(db_dir_path):
	#The 'Autometa databases' directory is empty
	print('AutoMeta databases directory empty, populating with appropriate databases.\nThis may take some time...')
	update_dbs(db_dir_path)

if not (os.path.isfile(db_dir_path + '/nr.dmnd') or os.path.isfile(db_dir_path + '/nr.dmnd.md5') or os.path.isfile(db_dir_path + '/nr.gz.md5')):
	print('NR database not found, downloading and building DIAMOND database.\nThis may take some time...')
	update_dbs(db_dir_path, 'nr')

if not (os.path.isfile(db_dir_path + '/prot.accession2taxid') or os.path.isfile(db_dir_path + '/prot.accession2taxid.gz.md5')):
	print('acc2taxid files not found, downloading.\nThis may take some time...')
	update_dbs(db_dir_path, 'acc2taxid')

if not (os.path.isfile(db_dir_path + '/names.dmp') or os.path.isfile(db_dir_path + '/nodes.dmp') or os.path.isfile(db_dir_path + '/taxdump.tar.gz.md5')):
	print('Taxdump files not found, downloading.\nThis may take some time...')
	update_dbs(db_dir_path, 'taxdump')

names_dmp_path = db_dir_path + '/names.dmp'
nodes_dmp_path = db_dir_path + '/nodes.dmp'
accession2taxid_path = db_dir_path + '/prot.accession2taxid'
diamond_db_path = db_dir_path + '/nr.dmnd'
current_taxdump_md5 = db_dir_path + '/taxdump.tar.gz.md5'
current_acc2taxid_md5 = db_dir_path + '/prot.accession2taxid.gz.md5'
current_nr_md5 = db_dir_path + '/nr.gz.md5'

if usr_prot_path:
	usr_prot_path = os.path.abspath(usr_prot_path)
	if os.path.isdir(usr_prot_path):
		print('You have provided a directory {}. \
		--user_prot_db requires a file path.'.format(usr_prot_path))
		exit(1)
	elif not os.path.isfile(usr_prot_path):
		print('{} is not a file.'.format(usr_prot_path))
		exit(1)
	else:
		diamond_db_path = usr_prot_path

if args['update']:
	print("Checking database directory for updates")
	update_dbs(db_dir_path, 'all')

if not os.path.isfile(prodigal_output + ".faa"):
	print "Prodigal output not found. Running prodigal..."
	#Check for file and if it doesn't exist run make_marker_table
	length_trim(fasta_path, length_cutoff)
	run_prodigal(filtered_assembly)

if not os.path.isfile(prodigal_output + ".daa"):
	print "Could not find {}. Running diamond blast... ".format(prodigal_output + ".daa")
	diamond_output = run_diamond(prodigal_output, diamond_db_path, num_processors, prodigal_daa)
elif os.stat(prodigal_output + ".daa").st_size == 0:
	print "{} file is empty. Re-running diamond blast...".format(prodigal_output + ".daa")
	diamond_output = run_diamond(prodigal_output, diamond_db_path, num_processors, prodigal_daa)
else:
	diamond_output = prodigal_output + ".tab"

if not os.path.isfile(prodigal_output + ".lca"):
	print "Could not find {}. Running lca...".format(prodigal_output + ".lca")
	blast2lca_output = run_blast2lca(diamond_output,db_dir_path)
elif os.stat(prodigal_output + ".lca").st_size == 0:
	print "{} file is empty. Re-running lca...".format(prodigal_output + ".lca")
	blast2lca_output = run_blast2lca(diamond_output,db_dir_path)
else:
	blast2lca_output = prodigal_output + ".lca"

print "Running add_contig_taxonomy.py... "
taxonomy_table = run_taxonomy(pipeline_path, filtered_assembly, blast2lca_output, db_dir_path, cov_table)

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
if not single_genome_mode:
	for kingdom in categorized_seq_objects:
		seq_list = categorized_seq_objects[kingdom]
		output_path = output_dir + '/' + kingdom + '.fasta'
		SeqIO.write(seq_list, output_path, 'fasta')

print "Done!"
