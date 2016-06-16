#!/usr/bin/env python

import argparse
import os
import subprocess 
import getpass
import time 
import multiprocessing
import pdb

#argument parser
parser = argparse.ArgumentParser(description="Script to run the autometa pipeline. \
	The script expects autometa repo to be somewhere in the user's home directory.")
parser.add_argument('-a','--assembly', help='assembly.fasta', required=True)
parser.add_argument('-p','--processors', help='assembly.fasta', default=1)
parser.add_argument('-l','--length_cutoff', help='Contig length cutoff to consider for binning.\
 Default is 10,000 bp.', default=10000, type = int)
args = vars(parser.parse_args())

length_cutoff = args['length_cutoff']
fasta_assembly = args['assembly']
processors = args['processors']
#kmer = args['kmer']

#def is_fasta(fasta):
#def process_assembly_name(fasta):
#check for output - so as not to run again

def length_trim(fasta,length_cutoff):
	#will need to update path of this perl script
	outfile_name = str(args['assembly'].split("/")[-1].split(".")[0]) + "_filtered.fasta"
	subprocess.call("{}fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff,outfile_name), shell = True)
	return outfile_name

def make_contig_table(fasta):
	#looks like this script is assuming contigs from a spades assembly
	output_table_name = str(fasta).split('.')[0] + ".tab"
	subprocess.call("{}make_contig_table.py {} {}".format(pipeline_path,fasta,output_table_name), shell = True)
	return output_table_name

def make_marker_table(fasta):
	hmm_marker_path = autometa_path + "/single-copy_markers/Bacteria_single_copy.hmm"
	hmm_cutoffs_path = autometa_path + "/single-copy_markers/Bacteria_single_copy_cutoffs.txt"
	#need to add processors to this script
	output_marker_table = fasta.split('.')[0] + "_marker.tab"
	if os.path.isfile(output_marker_table):
		print "{} file already exists!".format(output_marker_table)
		print "Continuing to next step..."
	else:
		print "Making the marker table with prodigal and hmmscan. This could take a while..."
		subprocess.call("hmmpress -f {}".format(hmm_marker_path), shell=True)
		subprocess.call("{}make_marker_table.py -a {} -m {} -c {} -o {} -p {}".\
			format(pipeline_path,fasta, hmm_marker_path, hmm_cutoffs_path,output_marker_table,args['processors']), \
			shell = True,stdout=FNULL, stderr=subprocess.STDOUT)
	return output_marker_table

def run_VizBin(fasta,marker_table):
	if os.path.isfile("contig_vizbin.tab"):
		print "VizBin output already exists!".format(marker_table)
		print "Continuing to next step..."
		return None
	else:
		subprocess.call("java -jar {}VizBin-dist.jar -i {} -o points.txt".format(autometa_path + "/VizBin/dist/",\
		fasta), shell = True)
		tmp_path = subprocess.check_output("ls /tmp/map* -dlt | grep {} | head -n1".format(username), shell=True).rstrip("\n").split()[-1]
		return tmp_path


def process_and_clean_VizBin(tmp_path,output_table_name):
	if tmp_path != None:
		subprocess.call("{}vizbin_process.pl {}/filteredSequences.fa points.txt > vizbin_table.txt".format(pipeline_path,tmp_path), shell=True)
		subprocess.call("tail -n +2 vizbin_table.txt | sort -k 1,1 > vizbin_table_sort.txt", shell=True)
		subprocess.call("tail -n +2 {} | sort -k 1,1 > contig_table_sort.txt".format(output_table_name), shell=True)
		subprocess.call("head -n 1 vizbin_table.txt > vizbin_header.txt", shell=True)
		subprocess.call("head -n 1 {} > contig_header.txt".format(output_table_name), shell=True)
		subprocess.call("join contig_header.txt vizbin_header.txt | sed $'s/ /\t/g' > joined_header.txt", shell=True)
		subprocess.call("touch contig_vizbin.tab; cat joined_header.txt >> contig_vizbin.tab; join contig_table_sort.txt vizbin_table_sort.txt |\
		 cat >> contig_vizbin.tab; sed $'s/ /\t/g' contig_vizbin.tab", shell=True,stdout=FNULL, stderr=subprocess.STDOUT)
		#Delete most recent /tmp/map* directory if it's the same user 
		subprocess.call("rm -rf {}".format(tmp_path), shell = True)
		#need more elegant way to clean up
		subprocess.call("ls -t *txt | head -n 7 | xargs -L1 rm".format(tmp_path), shell = True)

def install_VizBin_executable(autometa_path,home_dir):
	#install config into home directory 
	subprocess.call("cp -R {} {}".format(autometa_path + "/VizBin/.vizbin", home_dir), shell = True)
		#change config path
	subprocess.call("sed -i {}.vizbin/config 's?/home/user/'{}'?g'".format(home_dir,home_dir), shell = True)

def bin_assess_and_pick_cluster(marker_tab, vizbin_output_path):
	#Need to check for and install "dbscan" and "docopt" dependency from the command line (with CRAN mirror 27 [USA: MI])
	subprocess.call("Rscript {}dbscan_batch.R {} 0.3 5".format(pipeline_path, vizbin_output_path), shell = True,stdout=FNULL, stderr=subprocess.STDOUT)
	print "Running dbscan..."
	subprocess.call("{}assess_clustering.py -s {} -d *.tab_eps* -o assess_clustering_output".format(pipeline_path,marker_tab, vizbin_output_path), shell = True,stdout=FNULL, stderr=subprocess.STDOUT)
	print("The best cluster is:")
	subprocess.call("{}pick_best_clustering.py -i assess_clustering_output".format(pipeline_path), shell = True)
	best_cluster_tab = subprocess.check_output("{}pick_best_clustering.py -i assess_clustering_output".format(pipeline_path), shell = True)
	subprocess.call("mkdir -p eps_test_dir", shell = True)
	subprocess.call("mv *.tab_eps* eps_test_dir", shell = True)
	return best_cluster_tab.rstrip("\n")

def extract_best_clusters(fasta,best_cluster_tab,marker_tab_path):
	#hmm_marker_path = autometa_path + "/single-copy_markers/Bacteria_single_copy.hmm"
	#hmm_cutoffs_path = autometa_path + "/single-copy_markers/Bacteria_single_copy_cutoffs.txt"
	subprocess.call("mkdir -p best_cluster_output_dir", shell = True)
	#subprocess.call("{}cluster_separate_and_analyze.pl --fasta {} --table {} --outputdir best_cluster_output_dir --hmmdb {} --cutoffs {}\
		#".format(pipeline_path,fasta,best_cluster_tab,hmm_marker_path,hmm_cutoffs_path), shell = True)
		#use cluster_completeness.py instead
	subprocess.call("{}cluster_completeness.py -f {} -d eps_test_dir/{} -c 'db.cluster' -o best_cluster_output_dir -m {}\
		".format(pipeline_path,fasta,best_cluster_tab,marker_tab_path), shell = True)

start_time = time.time()
FNULL = open(os.devnull, 'w')

#check if fasta in path
if not os.path.isfile(args['assembly']):
	print "Could not find {}...".format(args['assembly'])
	exit()
#else:
	#check if fasta

#Check user CPUs
user_CPU_number = multiprocessing.cpu_count()

username = getpass.getuser()
home = os.path.expanduser("~") + "/"
autometa_path = subprocess.check_output('find ~ -name "autometa"', shell=True).rstrip("\n")
pipeline_path = autometa_path + "/pipeline/"
#Alternatively, the user could set this as an env variable


#run length trim and store output name
filtered_assembly = length_trim(args['assembly'],args['length_cutoff'])
contig_table = make_contig_table(filtered_assembly)
marker_tab_path = make_marker_table(filtered_assembly)
vizbin_output_path = "contig_vizbin.tab"

#install_VizBin_executable(autometa_path,home)

process_and_clean_VizBin(run_VizBin(filtered_assembly,marker_tab_path),contig_table)
#extract_best_clusters("scaffolds_over3k_over10k.fasta",bin_assess_and_pick_cluster("scaffolds_over3k_marker.tab", "contig_vizbin.tab"))
best_cluster_tab = bin_assess_and_pick_cluster(marker_tab_path, vizbin_output_path)
extract_best_clusters(filtered_assembly,best_cluster_tab,marker_tab_path)

elapsed_time = (time.time() - start_time)

print "Done!"
print "Elapsed time is {} seconds".format(round(elapsed_time,2))
FNULL.close()

