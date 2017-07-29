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

def length_trim(fasta,length_cutoff):
	#will need to update path of this perl script
	outfile_name = str(args['a'].split("/")[-1].split(".")[0]) + "_filtered.fasta"
	logger.info("{}/fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff,outfile_name))
	subprocess.call("{}/fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff,outfile_name), shell = True)
	return outfile_name

def make_contig_table(fasta):
	#looks like this script is assuming contigs from a spades assembly
	output_table_name = str(fasta).split('.')[0] + ".tab"
	logger.info("{}/make_contig_table.py {} {}".format(pipeline_path,fasta,output_table_name))
	subprocess.call("{}/make_contig_table.py {} {}".format(pipeline_path,fasta,output_table_name), shell = True)
	return output_table_name

def make_marker_table(fasta):
	if kingdom == 'bacteria':
		hmm_marker_path = autometa_path + "/single-copy_markers/Bacteria_single_copy.hmm"
		hmm_cutoffs_path = autometa_path + "/single-copy_markers/Bacteria_single_copy_cutoffs.txt"
	elif kingdom == 'archaea':
		hmm_marker_path = autometa_path + '/single-copy_markers/Archaea_single_copy.hmm'
		hmm_cutoffs_path = autometa_path + '/single-copy_markers/Archaea_single_copy_cutoffs.txt'

	#need to add processors to this script
	output_marker_table = fasta.split('.')[0] + "_marker.tab"
	if os.path.isfile(output_marker_table):
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
			format(pipeline_path,fasta, hmm_marker_path, hmm_cutoffs_path,output_marker_table, args['p']))
		subprocess.call("{}/make_marker_table.py -a {} -m {} -c {} -o {} -p {}".\
			format(pipeline_path,fasta, hmm_marker_path, hmm_cutoffs_path,output_marker_table, args['p']), \
			shell = True, stdout=FNULL, stderr=subprocess.STDOUT)
	return output_marker_table

def bin_assess_and_pick_cluster(pipeline_path,marker_tab, filtered_assembly, contig_table):
	if taxonomy_table_path:
		logger.info("{}/recursive_dbscan.py -m {} -d {} -f {} -o ./ -c {} -p {} -t {}".format(pipeline_path,marker_tab,kingdom,filtered_assembly, contig_table, processors, taxonomy_table_path))
		subprocess.call("{}/recursive_dbscan.py -m {} -d {} -f {} -o ./ -c {} -p {} -t {}".format(pipeline_path,marker_tab,kingdom,filtered_assembly, contig_table, processors, taxonomy_table_path), shell=True)
	else:
		logger.info("{}/recursive_dbscan.py -m {} -d {} -f {} -o ./ -c {} -p {}".format(pipeline_path,marker_tab,kingdom,filtered_assembly, contig_table, processors))
		subprocess.call("{}/recursive_dbscan.py -m {} -d {} -f {} -o ./ -c {} -p {}".format(pipeline_path,marker_tab,kingdom,filtered_assembly, contig_table, processors), shell=True)

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
parser.add_argument('-a', metavar='assembly', help='assembly.fasta', required=True)
parser.add_argument('-p', metavar='processors', help='Processors to use', default=1)
parser.add_argument('-l', metavar='cutoff length', help='Contig length cutoff to consider for binning.\
 Default is 10,000 bp.', default=10000, type = int)
parser.add_argument('-c', metavar='cluster completeness output', help='Best cluster output limited by completeness', default=20)
parser.add_argument('-k', metavar='kingdom', help='Kingdom to consider (archaea|bacteria)',\
choices=['bacteria','archaea'], default = 'bacteria')
parser.add_argument('-t', metavar='taxonomy table', help='Output of make_taxonomy_table.py')
args = vars(parser.parse_args())

length_cutoff = args['l']
fasta_assembly = os.path.basename(args['a'])
processors = args['p']
cluster_completeness = args['c']
kingdom = args['k'].lower()
taxonomy_table_path = args['t']

#check if fasta in path
if not os.path.isfile(args['a']):
	print "Could not find {}...".format(args['a'])
	logger.debug('Could not find {}...'.format(args['a']))
	exit()

#what input variables were and when you ran it (report fill path based on argparse)
logger.info('Input: -a {} -p {} -l {} -c {}'.format(args['a'], processors, length_cutoff, cluster_completeness))

start_time = time.time()
FNULL = open(os.devnull, 'w')

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

#run length trim and store output name
filtered_assembly = length_trim(args['a'], args['l'])
contig_table = make_contig_table(filtered_assembly)
marker_tab_path = make_marker_table(filtered_assembly)
vizbin_output_path = "contig_vizbin.tab"

bin_assess_and_pick_cluster(pipeline_path, marker_tab_path, filtered_assembly, contig_table)

elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))

print "Done!"
print "Elapsed time is {} (HH:MM:SS)".format(elapsed_time)
logger.info('Done!')
logger.info('Elapsed time is {} (HH:MM:SS)'.format(elapsed_time))
FNULL.close()
