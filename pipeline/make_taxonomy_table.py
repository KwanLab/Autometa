#!/usr/bin/env python

import subprocess
import argparse
import os
import pandas as pd
from Bio import SeqIO
import sys

def update_dbs(database_path, db='all'):
    """Updates databases for AutoMeta usage"""
    taxdump_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    taxdump_md5 = taxdump_url+".md5"
    accession2taxid_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
    accession2taxid_md5 = accession2taxid_url+".md5"
    nr_diamond_db_url = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
    nr_diamond_db_md5 = nr_diamond_db_url+".md5"
    #Downloading files for db population
    if db == 'all':
        os.system('wget -q %s -O %s/taxdump.tar.gz' % (taxdump_url, database_path))
        os.system('wget -q %s -O %s/taxdump.tar.gz.md5' % (taxdump_md5, database_path))
        os.system('tar -xzf %s/taxdump.tar.gz -C %s names.dmp nodes.dmp' % (database_path, database_path))
        os.system('rm %s/taxdump.tar.gz' % database_path)
        os.system('wget -q %s -O %s/prot.accession2taxid.gz' % (accession2taxid_url, database_path))
        os.system('wget -q %s -O %s/prot.accession2taxid.gz.md5' % (accession2taxid_md5, database_path))
        print("Gunzipping prot.accession2taxid gzipped file\nThis may take some time...")
        os.system('gunzip -9vNf %s/prot.accession2taxid.gz > %s/prot.accession2taxid' % (database_path, database_path))
        print("Downloading and building nr.dmnd database\nThis may take some time...")
        os.system('wget -q %s -O %s/nr.gz' % (nr_diamond_db_url, database_path))
        os.system('wget -q %s -O %s/nr.gz.md5' % (nr_diamond_db_md5, database_path))
        subprocess.call("diamond makedb --in {} --db {}/nr".format(database_path+'/nr.gz', database_path), shell = True)
        os.system('rm %s/nr.gz' % database_path)
        print("\nnodes.dmp, names.dmp, prot.accession2taxid and nr.dmnd files updated in %s\nResuming..." % database_path)
    if db == 'nr':
        print("updating nr.dmnd")
        os.system('wget -q %s -O %s/nr.gz' % (nr_diamond_db_url, database_path))
        os.system('wget -q %s -O %s/nr.gz.md5' % (nr_diamond_db_md5, database_path))
        print("building nr.dmnd database, this may take some time")
        subprocess.call("diamond makedb --in {} --db {}/nr".format(database_path+'/nr.gz', database_path), shell = True)
        os.system('rm %s/nr.gz' % database_path)
        print("nr.dmnd updated")
    if db == 'acc2taxid':
        print("updating prot.accession2taxid")
        os.system('wget -q %s -O %s.gz' % (accession2taxid_url, accession2taxid_path))
        os.system('wget -q %s -O %s.gz.md5' % (accession2taxid_md5, accession2taxid_path))
        print("Gunzipping prot.accession2taxid gzipped file\nThis may take some time...")
        os.system('gunzip -9vNf %s.gz > %s' % (accession2taxid_path, accession2taxid_path))
        print("prot.accession2taxid updated")
    if db == 'taxdump':
        print("updating nodes.dmp and names.dmp")
        os.system('wget -q %s -O %s/taxdump.tar.gz' % (taxdump_url, database_path))
        os.system('wget -q %s -O %s/taxdump.tar.gz.md5' % (taxdump_md5, database_path))
        os.system('tar -xzf %s/taxdump.tar.gz -C %s names.dmp nodes.dmp' % (database_path, database_path))
        os.system('rm %s/taxdump.tar.gz' % database_path)
        print("nodes.dmp and names.dmp updated")

def check_db(current_nr_md5, current_acc2taxid_md5, current_taxdump_md5):
    taxdump_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
    taxdump_md5 = taxdump_url+".md5"
    accession2taxid_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
    accession2taxid_md5 = accession2taxid_url+".md5"
    nr_diamond_db_url = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
    nr_diamond_db_md5 = nr_diamond_db_url+".md5"
    nr_status = os.system('wget -q -O- {} | md5sum | sed \'s:-$:{}:\' | md5sum -c --status'.format(nr_diamond_db_md5, current_nr_md5))
    acc2taxid_status = os.system('wget -q -O- {} | md5sum | sed \'s:-$:{}:\' | md5sum -c --status'.format(accession2taxid_md5, current_acc2taxid_md5))
    taxdump_status = os.system('wget -q -O- {} | md5sum | sed \'s:-$:{}:\' | md5sum -c --status'.format(taxdump_md5, current_taxdump_md5))
    return(nr_status, acc2taxid_status, taxdump_status)

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

def run_diamond(prodigal_output, diamond_db_path, num_processors, prodigal_daa):
    view_output = prodigal_output + ".tab"
    current_dir = os.getcwd()
    tmp_dir_path = current_dir + '/tmp'
    if not os.path.isdir(tmp_dir_path):
        os.makedirs(tmp_dir_path) # This will give an error if the path exists but is a file instead of a dir
    subprocess.call("diamond blastp --query {}.faa --db {} --evalue 1e-5 --max-target-seqs 200 -p {} --daa {} -t {}".format(prodigal_output, diamond_db_path, num_processors, prodigal_daa,tmp_dir_path), shell = True)
    subprocess.call("diamond view -a {} -f tab -o {}".format(prodigal_daa, view_output), shell = True)
    return view_output

#blast2lca using accession numbers#
def run_blast2lca(input_file, taxdump_path):
	output = input_file.rstrip(".tab") + ".lca"
	if os.path.isfile(output):
		print "{} file already exists!".format(output)
		print "Continuing to next step..."
	else:
		subprocess.call("{}/lca.py database_directory {} {} > {}".format(pipeline_path, db_dir_path, input_file, output), shell = True)
	return output

def run_taxonomy(pipeline_path, assembly_path, tax_table_path, db_dir_path): #Have to update this
	initial_table_path = assembly_path + '.tab'
	subprocess.call("{}/make_contig_table.py {} {}".format(pipeline_path, assembly_path, initial_table_path), shell = True)
	subprocess.call("{}/add_contig_taxonomy.py {} {} {} taxonomy.tab".format(pipeline_path, initial_table_path, tax_table_path, db_dir_path), shell = True)
	#contig_table_path, tax_table_path, db_dir_path, output_file_path
	return 'taxonomy.tab'

#argument parser
parser = argparse.ArgumentParser(description="Script to generate the contig taxonomy table.", epilog="Output will be directed to recursive_dbscan.py")
parser.add_argument('-a', metavar='assembly', help='assembly.fasta', required=True)
parser.add_argument('-p', metavar='processors', help='Processors to use. If not specified, will infer', default=1)
parser.add_argument('-db', metavar='databases directory', help='Path to directory with taxdump, protein accessions and diamond (NR) protein files', required=False, default='autometa_databases_directory')
parser.add_argument('-l', metavar='cutoff length', help='Contig length cutoff to consider for binning.\
 Default is 10,000 bp.', default=10000, type = int)
parser.add_argument("-update", required=False, action='store_true',\
 help='Checks/Adds/Updates: nodes.dmp, names.dmp, accession2taxid, nr.dmnd files within specified directory. Default will update to Autometa databases directory')

args = vars(parser.parse_args())

db_dir_path = args['db'].rstrip('/')
num_processors = args['p']
length_cutoff = args['l']
fasta_path = args['a']
fasta_assembly_prefix = os.path.splitext(os.path.basename(args['a']))[0]
prodigal_output = fasta_assembly_prefix + "_filtered.orfs"
prodigal_daa = prodigal_output + ".daa"
pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)
#add_contig_path = pipeline_path
filtered_assembly = fasta_assembly_prefix + "_filtered.fasta"

"""Sets current db files within dir specified by user or defaults to AutoMeta db dir"""
if db_dir_path == 'autometa_databases_directory':
    #Default settings
    db_dir_path = autometa_path + '/databases'
    if not os.path.isdir(db_dir_path):
        #Verify the 'Autometa databases' directory exists
        print('No databases directory found, creating and populating AutoMeta databases directory\nThis may take some time...')
        os.system('mkdir {}'.format(db_dir_path))
        update_dbs(db_dir_path)
    elif not os.listdir(db_dir_path):
        #The 'Autometa databases' directory is empty
        print('AutoMeta databases directory empty, populating with appropriate databases.\nThis may take some time...')
        update_dbs(db_dir_path)

names_dmp_path = db_dir_path + '/names.dmp'
nodes_dmp_path = db_dir_path + '/nodes.dmp'
accession2taxid_path = db_dir_path + '/prot.accession2taxid'
diamond_db_path = db_dir_path + '/nr.dmnd'
current_taxdump_md5 = db_dir_path + '/taxdump.tar.gz.md5'
current_acc2taxid_md5 = db_dir_path + '/prot.accession2taxid.gz.md5'
current_nr_md5 = db_dir_path + '/nr.gz.md5'

if args['update']:
    print("Checking database directory for updates")
    nr_status, acc2taxid_status, taxdump_status = check_db(current_nr_md5, current_acc2taxid_md5, current_taxdump_md5)
    if nr_status != 0:
        update_dbs(db_dir_path, db='nr')
    else:
        print("nr.dmnd already up to date")
    if acc2taxid_status != 0:
        update_dbs(db_dir_path, db='acc2taxid')
    else:
        print("prot.accession2taxid already up to date")
    if taxdump_status != 0:
        update_dbs(db_dir_path, db='taxdump')
    else:
        print("nodes.dmp and names.dmp already up to date")
    print('Resuming make_taxonomy_table.py...')

if not os.path.isfile(prodigal_output + ".faa"):
	print "Prodigal output not found. Running prodigal..."
	#Check for file and if it doesn't exist run make_marker_table
	length_trim(fasta_path, fasta_assembly_prefix, length_cutoff)
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
taxonomy_table = run_taxonomy(pipeline_path, filtered_assembly, blast2lca_output, db_dir_path)

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
