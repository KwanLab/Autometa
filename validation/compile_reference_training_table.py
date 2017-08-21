#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Script to tabulate reference training table (training based on marker contigs).')
parser.add_argument('-t','--contig_tab', help='Name of master contig table file', required=True)
parser.add_argument('-o','--out_tab', help='Name of output table file', required=True)
args = vars(parser.parse_args())

#Expects to have reference genomne in master table
master_table = pd.read_csv(args['contig_tab'], sep = "\t")

reference_training_list = []
for count,row in master_table.iterrows():
    contig = row['contig']
    reference_genome = row['reference_genome']
    num_markers = row['num_single_copies']
    if reference_genome != "misassembled" and num_markers >= 1:
        reference_training_list.append(reference_genome)
    else:
        reference_training_list.append("unclustered")

master_table['reference_training'] = reference_training_list

master_table.to_csv(args['out_tab'], sep='\t')
