#!/usr/bin/env python2.7

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
import sys
import urllib2
import subprocess
import os

import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO

def run_command(command_string, stdout_path = None):
	# Function that checks if a command ran properly
	# If it didn't, then print an error message then quit
	print('make_taxonomy_table.py, run_command: ' + command_string)
	if stdout_path:
		f = open(stdout_path, 'w')
		exit_code = subprocess.call(command_string, stdout=f, shell=True)
		f.close()
	else:
		exit_code = subprocess.call(command_string, shell=True)

	if exit_code != 0:
		print('make_taxonomy_table.py: Error, the command:')
		print(command_string)
		print('failed, with exit code ' + str(exit_code))
		exit(1)

def run_command_return(command_string, stdout_path = None):
	# Function that checks if a command ran properly.
	# If it didn't, then print an error message then quit
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
	filename = os.path.basename(file_url)
	md5name = os.path.basename(md5_url)

	while True:
		run_command('wget {} -O {}'.format(file_url, destination_dir + '/' + filename))
		run_command('wget {} -O {}'.format(md5_url, destination_dir + '/' + md5name))

		downloaded_md5 = subprocess.check_output(['md5sum', destination_dir + '/' + filename]).split(' ')[0]

		with open(destination_dir + '/' + md5name, 'r') as check_md5_file:
			check_md5 = check_md5_file.readline().split(' ')[0]

		if downloaded_md5 == check_md5:
			print('md5 checksum successful. Continuing...')
			break
		else:
			print('md5 checksum unsuccessful. Retrying...')

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
				print("md5 is not current. Updating nr.dmnd")
				download_file(database_path, nr_db_url, nr_db_md5_url)
		else:
			print("updating nr.dmnd")
			download_file(database_path, nr_db_url, nr_db_md5_url)

		# Now we make the diamond database
		print("building nr.dmnd database, this may take some time")
		returnCode = subprocess.call("diamond makedb --in {} --db {}/nr -p {}"\
			.format(database_path+'/nr.gz', database_path, num_processors), shell = True)
		if returnCode == 0: # i.e. job was successful
		#Make an md5 file to signal that we have built the database successfully
			os.remove('{}/nr.gz'.format(database_path))
			print("nr.dmnd updated")
	if db == 'all' or db == 'acc2taxid':
		# Download prot.accession2taxid.gz only if the version we have is not current
		if os.path.isfile(database_path + '/prot.accession2taxid.gz.md5'):
			if not md5IsCurrent(database_path + '/prot.accession2taxid.gz.md5', accession2taxid_md5_url):
				print("md5 is not current. Updating prot.accession2taxid")
				download_file(database_path, accession2taxid_url, accession2taxid_md5_url)
		else:
			print("updating prot.accession2taxid")
			download_file(database_path, accession2taxid_url, accession2taxid_md5_url)

		if os.path.isfile(database_path + '/prot.accession2taxid.gz'):
			print("Gunzipping prot.accession2taxid gzipped file\nThis may take some time...")
			run_command('gunzip -9vNf {}'\
				.format(database_path+'/prot.accession2taxid.gz'), database_path+'/prot.accession2taxid')
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
			os.remove('{}/taxdump.tar.gz'.format(database_path))
			print("nodes.dmp and names.dmp updated")

def check_dbs(db_path):
	'''
	Determines what files need to be downloaded/updated depending on
	what database path was specified
	'''
	if db_path == autometa_path + '/databases':
		db_dict = {
			'nr': ['nr.dmnd','nr.dmnd.md5','nr.gz.md5'],
			'acc2taxid': ['prot.accession2taxid.gz.md5','prot.accession2taxid'],
			'taxdump': ['names.dmp','nodes.dmp','taxdump.tar.gz.md5']
			}
	else:
		db_dict = {
			'nr': ['nr.dmnd'],
			'acc2taxid': ['prot.accession2taxid'],
			'taxdump': ['names.dmp','nodes.dmp']
			}
	db_files = os.listdir(db_path)
	for db in db_dict:
		for fh in db_dict[db]:
			if fh not in db_files:
				print('{0} database not found, downloading/formatting.\n\
				This may take some time...'.format(db))
				update_dbs(db_path, db)

def length_trim(fasta_path,length_cutoff):
	input_fname, ext = os.path.splitext(os.path.basename(fasta_path))
	#Trim the length of fasta file
	outfile_path = output_dir + '/' + input_fname + ".filtered" + ext
	run_command("{}/fasta_length_trim.pl {} {} {}"\
	.format(pipeline_path, fasta_path, length_cutoff, outfile_path))
	return outfile_path

def run_prodigal(path_to_assembly):
	assembly_fname, _ = os.path.splitext(os.path.basename(path_to_assembly))
	output_path = output_dir + '/' + assembly_fname + '.orfs.faa'
	if os.path.isfile(output_path):
		print "{} file already exists!".format(output_path)
		print "Continuing to next step..."
	else:
		run_command('prodigal -i {} -a {}/{}.orfs.faa -p meta -m -o {}/{}.txt'\
		.format(path_to_assembly, output_dir, assembly_fname, output_dir, assembly_fname))

def run_diamond(prodigal_output, diamond_db_path, num_processors, prodigal_diamond):
	view_output = prodigal_output + ".blastp"
	tmp_dir_path = os.path.dirname(prodigal_output) + '/tmp'
	if not os.path.isdir(tmp_dir_path):
		os.makedirs(tmp_dir_path) # This will give an error if the path exists but is a file instead of a dir
	error = run_command_return("diamond blastp --query {0}.faa --db {1} \
	--evalue 1e-5 --max-target-seqs 200 -p {2} --outfmt 6 --out {3} -t {4}"\
	.format(prodigal_output, diamond_db_path, num_processors, prodigal_diamond, tmp_dir_path))

	# If there is an error, attempt to rebuild NR
	if error == 134 or error == str(134):
		print('Fatal: Not enough disk space for diamond alignment archive!')
		exit(1)
	if error:
		print('Error when performing diamond blastp:\n{}\nAttempting to correct by rebuilding nr...'.format(error))
		update_dbs(db_dir_path, 'nr')

		run_command("diamond blastp --query {0}.faa --db {1} --evalue 1e-5 \
		--max-target-seqs 200 -p {2} --outfmt 6 --out {3} -t {4}"\
		.format(prodigal_output, diamond_db_path, num_processors, prodigal_diamond, tmp_dir_path))

	run_command("diamond view -a {} -f tab -o {}".format(prodigal_daa, view_output))

	return view_output

#blast2lca using accession numbers#
def run_blast2lca(input_file, taxdump_path):
	fname = os.path.splitext(os.path.basename(input_file))[0] + ".lca"
	output = '/'.join([output_dir, fname])
	if os.path.isfile(output) and not os.stat(output).st_size == 0:
		print "{} file already exists! Continuing to next step...".format(output)
	else:
		run_command("{0}/lca.py database_directory {1} {2}"\
			.format(pipeline_path, db_dir_path, input_file))
	return output

def run_taxonomy(pipeline_path, assembly_path, tax_table_path, db_dir_path,
		coverage_table, bgcs_path=None, orfs_path=None): #Have to update this
	assembly_fname, _ = os.path.splitext(os.path.basename(assembly_path))
	initial_table_path = output_dir + '/' + assembly_fname + '.tab'

	# Only make the contig table if it doesn't already exist
	if not os.path.isfile(initial_table_path):
		if coverage_table:
			run_command("{}/make_contig_table.py -a {} -o {} -c {}"\
			.format(pipeline_path, assembly_path, initial_table_path, coverage_table))
		elif single_genome_mode:
			run_command("{}/make_contig_table.py -a {} -o {} -n"\
			.format(pipeline_path, assembly_path, initial_table_path))
		else:
			run_command("{}/make_contig_table.py -a {} -o {}"\
			.format(pipeline_path, assembly_path, initial_table_path))
	if bgcs_path:
		run_command("{}/mask_bgcs.py2.7 --bgc {} --orfs {} --lca {}"
		.format(pipeline_path, bgcs_path, orfs_path, tax_table_path))
		unmasked_fname = os.path.basename(tax_table_path).replace('.lca', '.unmasked.tsv')
		tax_table_path = '{}/{}'.format(output_dir, unmasked_fname)
	#two_files_generated: *.masked.tsv, *.unmasked.tsv
	run_command("{}/add_contig_taxonomy.py {} {} {} {}/taxonomy.tab"\
	.format(pipeline_path, initial_table_path, tax_table_path, db_dir_path, output_dir))

	return output_dir + '/taxonomy.tab'

pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)

#argument parser
parser = ArgumentParser(description="Script to generate the contig taxonomy table.",
	epilog="Output will be directed to recursive_dbscan_output.tab")
parser.add_argument('-a', '--assembly', metavar='<assembly.fasta>',
	help='Path to metagenomic assembly fasta', required=True)
parser.add_argument('-p', '--processors', metavar='<int>',
	help='Number of processors to use.', type=int, default=1)
parser.add_argument('-db', '--db_dir', metavar='<dir>',
	help='Path to directory with taxdump, protein accessions and diamond (NR) \
	protein files. If this path does not exist, will create and download files.',
	required=False, default=autometa_path + '/databases')
parser.add_argument('-udb', '--user_prot_db', metavar='<user_prot_db>',
	help='Replaces the default diamond database (nr.dmnd)', required=False)
parser.add_argument('-l', '--length_cutoff', metavar='<int>',
	help='Contig length cutoff to consider for binning in bp', default=10000, type=int)
parser.add_argument('-v', '--cov_table', metavar='<coverage.tab>',
	help="Path to coverage table made by calculate_read_coverage.py. If this is \
	not specified then coverage information will be extracted from contig names (SPAdes format)",
	required=False)
parser.add_argument('-o', '--output_dir', metavar='<dir>',
	help='Path to directory to store output files', default='.')
parser.add_argument('-bgc', '--bgcs_dir', metavar='<dir>',
	help='Path to directory of biosynthetic gene clusters. Masks BGCs')
parser.add_argument('-s', '--single_genome', help='Specifies single genome mode',
	action='store_true')
parser.add_argument('-u', '--update', required=False, action='store_true',
	help='Checks/Adds/Updates: nodes.dmp, names.dmp, accession2taxid, nr.dmnd files within specified directory.')

args = vars(parser.parse_args())

db_dir_path = os.path.abspath(args['db_dir'])
usr_prot_path = args['user_prot_db']
num_processors = args['processors']
length_cutoff = args['length_cutoff']
fasta_path = args['assembly']
cov_table = args['cov_table']
output_dir = args['output_dir']
single_genome_mode = args['single_genome']

bgcs_dir = args['bgcs_dir']
fasta_fname, _ = os.path.splitext(os.path.basename(fasta_path))
prodigal_output = '/'.join([output_dir, fasta_fname + ".filtered.orfs"])
prodigal_diamond = prodigal_output + ".blastp"

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
	print('No databases directory found, creating and populating AutoMeta databases directory\n\
	This may take some time...')
	os.mkdir(db_dir_path)
	update_dbs(db_dir_path)
elif not os.listdir(db_dir_path):
	#The 'Autometa databases' directory is empty
	print('AutoMeta databases directory empty, populating with appropriate databases.\n\
	This may take some time...')
	update_dbs(db_dir_path)
else:
	check_dbs(db_dir_path)

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

filtered_assembly = output_dir + '/' + fasta_fname + ".filtered.fasta"
if not os.path.isfile(filtered_assembly):
	filtered_assembly = length_trim(fasta_path, length_cutoff)

if not os.path.isfile(prodigal_output + ".faa"):
	print "Prodigal output not found. Running prodigal..."
	#Check for file and if it doesn't exist run make_marker_table
	run_prodigal(filtered_assembly)

if not os.path.isfile(prodigal_diamond):
	print "Could not find {}. Running diamond blast... ".format(prodigal_diamond)
	diamond_output = run_diamond(prodigal_output, diamond_db_path, num_processors, prodigal_diamond)
elif os.stat(prodigal_diamond).st_size == 0:
	print "{} file is empty. Re-running diamond blast...".format(prodigal_diamond)
	diamond_output = run_diamond(prodigal_output, diamond_db_path, num_processors, prodigal_diamond)
elif not os.path.isfile(prodigal_diamond):
	print "{0} not found. \nExiting..."\
 	.format(prodigal_diamond)
	exit(1)
else:
	diamond_output = prodigal_diamond

if not os.path.isfile(prodigal_output + ".lca"):
	print "Could not find {}. Running lca...".format(prodigal_output + ".lca")
	blast2lca_output = run_blast2lca(diamond_output,db_dir_path)
elif os.stat(prodigal_output + ".lca").st_size == 0:
	print "{} file is empty. Re-running lca...".format(prodigal_output + ".lca")
	blast2lca_output = run_blast2lca(diamond_output,db_dir_path)
else:
	blast2lca_output = prodigal_output + ".lca"

taxonomy_table = output_dir + '/taxonomy.tab'
if not os.path.isfile(taxonomy_table) or os.stat(taxonomy_table).st_size == 0:
	print "Running add_contig_taxonomy.py... "
	if bgcs_dir:
		taxonomy_table = run_taxonomy(pipeline_path=pipeline_path,
			assembly_path=filtered_assembly,
			tax_table_path=blast2lca_output,
			db_dir_path=db_dir_path,
			coverage_table=cov_table,
			bgcs_path=bgcs_dir,
			orfs_path=prodigal_output + '.faa')
	else:
		taxonomy_table = run_taxonomy(pipeline_path=pipeline_path,
			assembly_path=filtered_assembly,
			tax_table_path=blast2lca_output,
			db_dir_path=db_dir_path,
			coverage_table=cov_table)
else:
	print('taxonomy.tab exists... Splitting original contigs into kingdoms')

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
	if contig not in all_seq_records:
		#Using filtered assembly, taxonomy.tab contains contigs not filtered
		print('{0} below length filter, skipping.'.format(contig))
		continue
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
