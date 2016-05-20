#!/usr/bin/env python

import pandas as pd 
import argparse
import subprocess


#argument parser
parser = argparse.ArgumentParser(description='Script tabulate single copy markers \
	from a metagenome assembly. Dependencies: prodigal v2.6.2, hhmscan (hmmer 3.1b2)')
parser.add_argument('-a','--assembly', help='Input assembly file', required=True)
parser.add_argument('-c','--cutoffs', help='Bacterial single copy hmm cutoffs as defined by Rinke et al.', default="~/Bacteria_single_copy_cutoffs.txt")
parser.add_argument('-m','--hmm', help='Bacteria_single_copy_cutoffs.hmm', default="~/Bacteria_single_copy.hmm")
parser.add_argument('-o','--out', help='outfile.tab for normal mixute modeling in R', required=True)
args = vars(parser.parse_args())

assembly = args['assembly']

def get_contig_list(path_to_assembly):
    #Get list of spades contigs
    assembly_handle = open(path_to_assembly,"rU")
    contig_name_list = []
    for line in assembly_handle:
        if '>' in line:
            contig_name = line.rstrip("\n").split()[0][1:]
            contig_name_list.append(contig_name)
    assembly_handle.close()
    return contig_name_list

def run_prodigal(path_to_assembly):
	#When "shell = True", need to give one string, not a list
	subprocess.call(" ".join(['prodigal ','-i ' + path_to_assembly, '-a ' + path_to_assembly.split(".")[0] +\
	 '.orfs.faa','-p meta' '-m', '-o ' + path_to_assembly.split(".")[0] + '.txt']), shell = True)

def run_hhmscan(path_to_prodigal_output,hmmdb):
	subprocess.call("hmmscan --tblout {} {} {}".format(path_to_prodigal_output.split(".")[0] + ".hmm.tbl", hmmdb, path_to_prodigal_output),shell = True)

run_prodigal(assembly)
run_hhmscan(assembly.split(".")[0] + '.orfs.faa',args['hmm'])

hmm_table = pd.read_csv(assembly.split(".")[0] + ".hmm.tbl", sep='\s+', \
	usecols = [1,2,5], skiprows = 3, header = None, index_col=False, engine = 'python')
cutoffs_table = pd.read_csv(args['cutoffs'], sep = '\s', engine = 'python', header = None)

#Search for contigs/ORFs that contain single copy PFAM domains that pass cutoff
#for loop to search for PFAM domains in hmm table column 1:
contig_ORFs_that_pass_cutoffs = {}
for index,PFAM_cutoffs_id in enumerate(cutoffs_table[0]):
	for count,PFAM_hmm_scan_id in enumerate(hmm_table[1]):	
		#If the PFAMs domain match (cutoffs ID names format are PFAMXXXXX, not PFAMXXXXX.1)
		if str(PFAM_cutoffs_id) in str(PFAM_hmm_scan_id):
			#and hmm score value > cutoff value, 
			if float(hmm_table[5][count]) > float(cutoffs_table[1][index]):
				#make dictionary of lists for PFAMs, where PFAM is key and contigs populate the list
				if PFAM_hmm_scan_id not in contig_ORFs_that_pass_cutoffs:
					contig_ORFs_that_pass_cutoffs[PFAM_hmm_scan_id] = [] 
				else:
					contig_ORFs_that_pass_cutoffs[PFAM_hmm_scan_id].append(hmm_table[2][count])

#Make a dictionary of dictionaries with contigs that have single copy genes (list PFAMS, for final table count length of list),
#contig length, contig GC, contig len, passecd PFAM domains
#write out tab-delimited table

with open(args['out'], 'w') as outfile:
	outfile.write("contig_name" + '\t'+ "single_copy_PFAMs" + '\t' + "num_single_copies" + '\n')
	contig_dictionary = {}
	for count,contig in enumerate(get_contig_list(assembly)):
		contig_dictionary[contig] = {}
		contig_dictionary[contig]['single_copy_PFAMs'] = []
		contig_dictionary[contig]['num_single_copies'] = 0
		for PFAM_key,contigs in contig_ORFs_that_pass_cutoffs.items():
			for item in contigs:
				if str(contig) in item:
					contig_dictionary[contig]['single_copy_PFAMs'].append(PFAM_key)
		contig_dictionary[contig]['num_single_copies'] = len(contig_dictionary[contig]['single_copy_PFAMs'])

		if len(contig_dictionary[contig]['single_copy_PFAMs']) > 0:
			outfile.write(str(contig) + '\t' + str(",".join(contig_dictionary[contig]['single_copy_PFAMs'])) + '\t' + \
		 	  str(contig_dictionary[contig]['num_single_copies']) + '\n')
		else:
			outfile.write(str(contig) + '\t' + "NA" + '\t' + \
		 	  str(contig_dictionary[contig]['num_single_copies']) + '\n')			

print("\nDone!")