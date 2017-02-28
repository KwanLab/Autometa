#!/usr/bin/env python

import pandas as pd
import argparse

#1. Load table into pandas df
#2. Populate dictionary, key = contig, value = reference genome
#3. Order output based on names in assembly file (make list of contig names)
#4. Report number of contigs that don't have an alignment

#Note - check for multiple alignments

parser = argparse.ArgumentParser(description='Script to assign reference genome\
    to de novo contig using nucmer alignment produced by quast')
parser.add_argument('-a','--assembly', help='De novo assembly nucl fasta file', required=True)
parser.add_argument('-c','--coordinate_table', help='Nucmer *.coords* table', required=True)
parser.add_argument('-f','--filtered', help='Nucmer .coords table is filtered', default=True)
parser.add_argument('-o','--out_table', help='Output.tsv', default="contig_reference.tsv")
args = vars(parser.parse_args())

coords_table = args['coordinate_table']
path_to_asm = args['assembly']
filter_boolean = args['filtered']

def parse_filtered_nucmer_table(coords_table,filtered=True):
    #1. Load table into pandas df
    filtered_coordinate_table = pd.read_csv(coords_table, header=None,skiprows=2,sep="|")
    #2. Populate dictionary, key = contig, value = reference genome
    alignment_dict = {}
    for count,alignment_info in enumerate(filtered_coordinate_table[4]):
        #Looks like the reference genome name is truncated to a certain # of chars.
        if filtered==False:
            #print alignment_info.lstrip(" ").split("\t")
            reference_genome = alignment_info.lstrip(" ").split("\t")[0]
            contig = alignment_info.lstrip(" ").split("\t")[1]
        else:
            reference_genome = alignment_info.lstrip(" ").split()[0]
            contig = alignment_info.lstrip(" ").split()[1]

        print contig,reference_genome
        alignment_dict[contig] = reference_genome
    return alignment_dict

alignment_dict = parse_filtered_nucmer_table(coords_table,filter_boolean)

#3. Order output based on names in assembly file (make list of contig names)
contig_list = []
with open(path_to_asm, "r") as asm_file:
    for line in asm_file:
        if ">" in line:
            contig_list.append(line.lstrip(">").rstrip())

#4. Report number of contigs that don't have an alignment
unaligned_contigs = []
with open("out.tab", "w") as outfile:
    for contig in contig_list:
        reference_genome = "NA"
        if contig in alignment_dict:
            reference_genome = alignment_dict[contig]
        else:
            unaligned_contigs.append(contig)
        outfile.write(contig + "\t" + reference_genome + "\n")

print("{} of {} contigs were unaligned by nucmer...".format(len(unaligned_contigs),len(contig_list)))
