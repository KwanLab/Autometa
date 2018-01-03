#!/usr/bin/env python

import argparse
import os
import sys
import subprocess
import getpass
import time
import multiprocessing
import pdb
import logging

def cythonize_lca_functions():
	logger.info("{}/lca_functions.so not found, cythonizing lca_function.pyx for make_taxonomy_table.py".format(pipeline_path))
	subprocess.call("cd {}".format(pipeline_path), shell=True)
	subprocess.call("./setup_lca_functions.py build_ext --inplace", shell = True)
	subprocess.call("cd {}".format(output_dir), shell=True)

def run_make_taxonomy_tab(fasta, length_cutoff):
	"""Runs make_taxonomy_table.py and directs output to taxonomy.tab for run_autometa.py"""
	# Note we don't have to supply the cov_table here because earlier in this script we already run make_contig_table.py
	output_path = output_dir + '/taxonomy.tab'
	logger.info("{}/make_taxonomy_table.py  -a {} -db {} -p {} -l {}".\
		format(pipeline_path, fasta, db_dir_path, processors, length_cutoff))
	subprocess.call("{}/make_taxonomy_table.py -a {} -db {} -p {} -l {}".\
		format(pipeline_path, fasta, db_dir_path, processors, length_cutoff),\
		shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
	return output_path

def length_trim(fasta,length_cutoff):
	#will need to update path of this perl script
	outfile_name = os.path.basename(fasta).split(".")[0] + "_filtered.fasta"
	output_path = output_dir + '/' + outfile_name
	logger.info("{}/fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff,output_path))
	subprocess.call("{}/fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff,output_path), shell = True)
	return outfile_name

def make_contig_table(fasta, coverage_table):
	output_table_name = str(fasta).split('.')[0] + ".tab"
	output_path = output_dir + '/' + output_table_name
	if coverage_table:
		logger.info("{}/make_contig_table.py -a {} -c {} -o {}".format(pipeline_path,fasta,coverage_table,output_path))
		subprocess.call("{}/make_contig_table.py -a {} -c {} -o {}".format(pipeline_path,fasta,coverage_table,output_path), shell = True)
	else:
		logger.info("{}/make_contig_table.py -a {} -o {}".format(pipeline_path,fasta,output_path))
		subprocess.call("{}/make_contig_table.py -a {} -o {}".format(pipeline_path,fasta,output_path), shell = True)
	return output_path

def make_marker_table(fasta):
	if kingdom == 'bacteria':
		hmm_marker_path = autometa_path + "/single-copy_markers/Bacteria_single_copy.hmm"
		hmm_cutoffs_path = autometa_path + "/single-copy_markers/Bacteria_single_copy_cutoffs.txt"
	elif kingdom == 'archaea':
		hmm_marker_path = autometa_path + '/single-copy_markers/Archaea_single_copy.hmm'
		hmm_cutoffs_path = autometa_path + '/single-copy_markers/Archaea_single_copy_cutoffs.txt'

	#need to add processors to this script
	output_marker_table = fasta.split('.')[0] + "_marker.tab"
	output_path = output_dir + '/' + output_marker_table
	if os.path.isfile(output_path):
		print "{} file already exists!".format(output_marker_table)
		print "Continuing to next step..."
		logger.info('{} file already exists!'.format(output_marker_table))
		logger.info('Continuing to next step...')
	else:
		print "Making the marker table with prodigal and hmmscan. This could take a while..."
		logger.info('Making the marker table with prodigal and hmmscan. This could take a while...')
		logger.info("hmmpress -f {}".format(hmm_marker_path))
		subprocess.call("hmmpress -f {}".format(hmm_marker_path), shell=True,stdout=FNULL, stderr=subprocess.STDOUT)
		logger.info("{}/make_marker_table.py -a {} -m {} -c {} -o {} -p {}".\
			format(pipeline_path,fasta, hmm_marker_path, hmm_cutoffs_path,output_path, args['p']))
		subprocess.call("{}/make_marker_table.py -a {} -m {} -c {} -o {} -p {}".\
			format(pipeline_path,fasta, hmm_marker_path, hmm_cutoffs_path,output_path, args['p']), \
			shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
	return output_path

def recursive_dbscan(input_table, filtered_assembly, domain):
	recursive_dbscan_output_path = output_dir + '/recursive_dbscan_output.tab'
	k_mer_file = output_dir + '/k-mer_matrix'
	logger.info("{}/recursive_dbscan.py -t {} -a {} -d {} -k {}".format(pipeline_path, input_table, filtered_assembly, output_dir, domain))
	subprocess.call("{}/recursive_dbscan.py -t {} -a {} -d {} -k {}".format(pipeline_path, input_table, filtered_assembly, output_dir, domain), shell=True)

	return recursive_dbscan_output_path, k_mer_file

def combine_tables(table1_path, table2_path):
	comb_table_path = output_dir + '/combined_contig_info.tab'
	# Note: in this sub we assume that the tables have column 1 in common
	# Store lines of table 2, keyed by the value of the first column

	# First make data structure from table 2
	table2_lines = dict()
	with open(table2_path) as table2:
		for i, line in enumerate(table2):
			line_list = line.rstrip().split('\t')
			contig = line_list.pop(0)
			if i == 0:
				table_2_header = '\t'.join(line_list)
			else:
				table2_lines[contig] = '\t'.join(line_list)

	comb_table = open(comb_table_path, 'w')
	with open(table1_path) as table1:
		for i, line in enumerate(table1):
			line_list = line.rstrip().split('\t')
			contig = line_list.pop(0)
			if i == 0:
				new_header = line.rstrip() + '\t' + table_2_header + '\n'
				comb_table.write(new_header)
			else:
				# We have to check whether the line exists in table 2
				if contig in table2_lines:
					new_line = line.rstrip() + '\t' + table2_lines[contig] + '\n'
					comb_table.write(new_line)

	comb_table.close()

	return comb_table_path

def ML_recruitment(input_table, matrix):
	ML_recruitment_output_path = output_dir + '/ML_recruitment_output.tab'
	logger.info("{}/ML_recruitment.py -t {} -p {} -r -m {} -o {}".format(pipeline_path, input_table, processors, matrix, ML_recruitment_output_path))
	subprocess.call("{}/ML_recruitment.py -t {} -p {} -r -m {} -o {}".format(pipeline_path, input_table, processors, matrix, ML_recruitment_output_path), shell=True)

	return ML_recruitment_output_path

#logger
logger = logging.getLogger('run_autometa.py')
hdlr = logging.FileHandler('run_autometa.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.DEBUG)

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
taxdump_dir = args['db_dir']
cov_table = args['cov_table']

#Check user CPUs
user_CPU_number = multiprocessing.cpu_count()

pipeline_path = sys.path[0]
pathList = pipeline_path.split('/')
pathList.pop()
autometa_path = '/'.join(pathList)

# Output current branch and commit
branch_command = "git -C " + autometa_path + " branch | grep \* | sed 's/^..//'"
branch = subprocess.Popen(branch_command, shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()

commit_command = 'git -C ' + autometa_path + ' rev-parse --short HEAD'
commit = subprocess.Popen(commit_command, shell=True, stdout=subprocess.PIPE).communicate()[0].rstrip()

logger.info('Currently running branch ' + branch + ', commit ' + commit)

#check if appropriate databases specified for make taxonomy table
if args['maketaxtable'] and not args['db']:
	print("Must specify databases directory (-db)")
	exit()

#check if fasta in path
if not os.path.isfile(args['a']):
	print "Could not find {}...".format(args['a'])
	logger.debug('Could not find {}...'.format(args['a']))
	exit()

#what input variables were and when you ran it (report fill path based on argparse)
logger.info('Input command: ' + ' '.join(sys.argv))

start_time = time.time()
FNULL = open(os.devnull, 'w')

#run length trim and store output name
filtered_assembly = length_trim(fasta_assembly, length_cutoff)
if cov_table:
	contig_table = make_contig_table(filtered_assembly, cov_table)
else:
	contig_table = make_contig_table(filtered_assembly, cov_table)
marker_tab_path = make_marker_table(filtered_assembly)

# Make combined table
if taxonomy_table_path and not args['maketaxtable']:
	combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
elif taxonomy_table_path and args['maketaxtable']:
	if not os.path.isfile(args['t']):
		print "Could not find {}, running make_taxonomy_table.py".format(args['t'])
		logger.debug('Could not find {}, running make_taxonomy_table.py'.format(args['t']))
		if not os.path.isfile(pipeline_path+"/lca_functions.so"):
			cythonize_lca_functions()
		taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
		combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
	elif os.path.isfile(args['t'] and os.stat(args['t']).st_size == 0):
		print "{} file is empty, running make_taxonomy_table.py".format(args['t'])
		logger.debug('{} file is empty, running make_taxonomy_table.py'.format(args['t']))
		if not os.path.isfile(pipeline_path+"/lca_functions.so"):
			cythonize_lca_functions()
		taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
		combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
	else:
		print "{} already exists, not performing make_taxonomy_table.py".format(args['t'])
		combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
elif not taxonomy_table_path and args['maketaxtable']:
	if not os.path.isfile(pipeline_path+"/lca_functions.so"):
		cythonize_lca_functions()
	taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
	combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
else:
	combined_table_path = combine_tables(contig_table, marker_tab_path)

recursive_dbscan_output, matrix_file = recursive_dbscan(combined_table_path, filtered_assembly, kingdom)

if do_ML_recruitment:
	ML_recruitment(recursive_dbscan_output, matrix_file)

elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))

print "Done!"
print "Elapsed time is {} (HH:MM:SS)".format(elapsed_time)
logger.info('Done!')
logger.info('Elapsed time is {} (HH:MM:SS)'.format(elapsed_time))
FNULL.close()
