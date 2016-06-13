#!/usr/bin/env python

import argparse
from os.path import expanduser
import subprocess 
import getpass

#argument parser
parser = argparse.ArgumentParser(description="Script to run the autometa pipeline. \
	The script expects autometa repo to be somewhere in the user's home directory.")
parser.add_argument('-a','--assembly', help='assembly.fasta', required=True)
parser.add_argument('-p','--processors', help='assembly.fasta', default=1)
parser.add_argument('-l','--length_cutoff', help='Contig length cutoff to consider for binning.\
 Default is 10,000 bp.', default=10000)
args = vars(parser.parse_args())

username = getpass.getuser()
home = expanduser("~") + "/"
autometa_path = subprocess.check_output('find ~ -name "autometa"', shell=True).rstrip("\n")
pipeline_path = autometa_path + "/pipeline/"
#Alternatively, the user could set this as an env variable

length_cutoff = args['length_cutoff']
fasta_assembly = args['assembly']
processors = args['processors']

processed_assembly_name = str(args['assembly'].split('.')[0]) + "_over{0}k.fasta".format(int(args['length_cutoff']/1000))
#def is_fasta(fasta):
#def process_assembly_name(fasta):
#check for output - so as not to run again

def length_trim(fasta,length_cutoff):
	#will need to update path of this perl script
	outfile_name = str(args['assembly'].split('.')[0]) + "_over{0}k.fasta".format(int(args['length_cutoff']/1000))
	subprocess.call("{}fasta_length_trim.pl {} {} {}".format(pipeline_path, fasta, length_cutoff,outfile_name), shell = True)

def make_contig_table(fasta):
	#looks like this script is assuming contigs from a spades assembly
	infile_name = str(args['assembly'].split('.')[0]) + "_over{0}k.fasta".format(int(args['length_cutoff']/1000))
	output_table_name = str(args['assembly'].split('.')[0]) + "_over{0}k.tab".format(int(args['length_cutoff']/1000))
	subprocess.call("{}make_contig_table.py {} {}".format(pipeline_path,infile_name, output_table_name), shell = True)

def make_marker_table(fasta):
	hmm_marker_path = autometa_path + "/single-copy_markers/Bacteria_single_copy.hmm"
	hmm_cutoffs_path = autometa_path + "/single-copy_markers/Bacteria_single_copy_cutoffs.txt"
	infile_name = str(args['assembly'].split('.')[0]) + "_over{0}k.fasta".format(int(args['length_cutoff']/1000))
	#need to add processors to this script
	subprocess.call("{}make_marker_table.py -a {} -m {} -c {}".format(pipeline_path,infile_name, hmm_marker_path, hmm_cutoffs_path), shell = True)

def run_and_process_VizBin(fasta):
	subprocess.call("java -jar {}VizBin-dist.jar -i {} -o points.txt".format(autometa_path + "/VizBin/dist/",\
		processed_assembly_name), shell = True)
	#process
	tmp_path = subprocess.check_output("ls /tmp/map* -dlt | grep {} | head -n1".format(username), shell=True).rstrip("\n").split()[-1]
	subprocess.call("vizbin_process.pl {}/filteredSequences.fa points.txt > vizbin_table".format(tmp_path), shell=True)
	subprocess.call("tail -n +2 vizbin_table | sort -k 1,1 > vizbin_table_sort", shell=True)
	contig_tab = str(args['assembly'].split('.')[0]) + "_over{0}k.tab".format(int(args['length_cutoff']/1000))
	subprocess.call("tail -n +2 {} | sort -k 1,1 > contig_table_sort".format(contig_tab), shell=True)
	subprocess.call("head -n 1 vizbin_table > vizbin_header", shell=True)
	subprocess.call("head -n 1 {} > contig_header".format(contig_tab), shell=True)
	subprocess.call("join contig_header vizbin_header | sed $'s/ /\t/g' > joined_header", shell=True)
	subprocess.call("touch contig_vizbin.tab; cat joined_header >> contig_vizbin.tab; join contig_table_sort vizbin_table_sort |\
	 cat >> contig_vizbin.tab; sed $'s/ /\t/g' contig_vizbin.tab", shell=True)
	#Delete most recent /tmp/map* directory if it's the same user 
	subprocess.call("rm -rf {}".format(tmp_path), shell = True)

def install_VizBin_executable(autometa_path,home_dir):
	subprocess.call("cp -R {} {}".format(autometa_path + "/VizBin/.vizbin", home_dir), shell = True)


length_trim(fasta_assembly,length_cutoff)
make_contig_table(fasta_assembly)
make_marker_table(fasta_assembly)
install_VizBin_executable(autometa_path,home)
run_VizBin(processed_assembly_name)

