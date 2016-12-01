#!/usr/bin/env python

import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Script to pull out fasta files from a contig table column")
parser.add_argument('-t','--contig_tab', help='Contig table', required=True)
parser.add_argument('-c','--cluster_column', help='Name of column for cluster', default='cluster')
parser.add_argument('-f','--fasta', help='Assembly FASTA file', required=True)
parser.add_argument('-o','--outdir', help='Path of directory for output', default='./')
args = vars(parser.parse_args())

contig_tab = pd.read_csv(args['contig_tab'], sep = "\t")
cluster_column = args['cluster_column']
fasta_path = args['fasta']
fasta_record_dict =  SeqIO.index(fasta_path, "fasta")

#Build cluster dict from contig table
cluster_dict = {}
for count,contig in enumerate(contig_tab['contig']):
    cluster = contig_tab[cluster_column][count]
    if cluster in cluster_dict:
        cluster_dict[cluster].append(contig)
    else:
        cluster_dict[cluster] = [contig]

#Make a list of seq records for each cluster and write out the fasta files
for cluster in cluster_dict:
    if cluster == 'unclustered':
        fasta_output_path = str(cluster) + '.fasta'
    else:
        fasta_output_path = 'cluster_' + str(cluster) + '.fasta'
    cluster_seq_record_list = cluster_dict[cluster]
    #Turn list of contig names into seq records
    seq_record_list = [fasta_record_dict[seq] for seq in cluster_seq_record_list]
    SeqIO.write(seq_record_list, fasta_output_path, 'fasta')
