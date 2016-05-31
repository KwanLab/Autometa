#!/usr/bin/env python

from __future__ import division
import pandas as pd 
import argparse


#argument parser
parser = argparse.ArgumentParser(description='Script that reads through \
	single copy contig tables can optionally output a histogram table. \
	Add the -v flag to see an estimate of # of genomes')
parser.add_argument('-t','--tab', help='Input marker .tab file', required=True)
parser.add_argument('-o','--out', help='optional outfile.tab for histogram', required=False)
parser.add_argument('-v','--verbose', help='prints out PFAM and corresponding frequency', required=False, action='store_true')
args = vars(parser.parse_args())

contig_table = pd.read_csv(args['tab'], sep = '\t', engine = 'python')

PFAM_dict = {}
for count,PFAMs in enumerate(contig_table['single_copy_PFAMs']):
	#could put a 'try' and 'except' loop here to check for propper formatting
	if str(PFAMs).startswith("PF"):
		for PFAM in str(PFAMs).split(','):
			#print PFAM
			#make dict of dict, where first level of key is PFAM, within each PFAM dict, there are contig size (len) and frequencies
			if PFAM in PFAM_dict:
				PFAM_dict[PFAM]['occurence'] += 1
			else:
				PFAM_dict[PFAM] = {}
				PFAM_dict[PFAM]['occurence'] = 1

if args['verbose'] == True:
	copy_average_sum = 0
	PFAM_frequency_list = []
	for PFAM,items in PFAM_dict.items():
		print PFAM,PFAM_dict[PFAM]['occurence']
		copy_average_sum += PFAM_dict[PFAM]['occurence']
		PFAM_frequency_list.append(PFAM_dict[PFAM]['occurence'])
	average_of_averages = (copy_average_sum / len(PFAM_dict.keys()))
	print("Average of single copy averages (estimate of # of genomes): %i") % (average_of_averages)
	print("Minimum is %i and maximum is %i") % (min(PFAM_frequency_list),max(PFAM_frequency_list))

if args['out'] != None:
	with open(str(args['out']), 'w') as outfile:
		outfile.write("PFAM" + '\t' + "Occurrences" + '\n')
		for PFAM,items in PFAM_dict.items():
			outfile.write(str(PFAM) + '\t' + str(PFAM_dict[PFAM]['occurence']) + '\n')