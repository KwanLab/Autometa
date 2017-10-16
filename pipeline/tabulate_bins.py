#!/usr/bin/env python

import glob
import argparse
from os import path

parser = argparse.ArgumentParser(description='Script to build contig-cluster table from fasta files.')
parser.add_argument('-i','--input_path', help='Path to directory with fasta files.', required=True)
parser.add_argument('-c','--cluster_name', help='Name of cluster column (Default: cluster)', default='cluster')
parser.add_argument('-x','--extension', help='Fasta file extension (Default: fasta)', default='fasta')
parser.add_argument('-o','--outfile', help='Output table path', default='contig_bins.tab')
args = vars(parser.parse_args())

dir_path = path.realpath(args['input_path'])
fasta_regex =  dir_path + "/*." + args['extension']

with open(args['outfile'],'w') as outfile:
    outfile.write("contig\t{}\n".format(args['cluster_name']))
    for fasta in glob.glob(fasta_regex):
        try:
            fasta_handle = open(fasta)
            for line in fasta_handle:
                if ">" in line:
                    contig_name = line.rstrip().lstrip(">")
                    bin_name = ".".join(path.basename(fasta).split(".")[:-1])
                    outfile.write(contig_name + "\t" + bin_name + "\n")
        except IOError as e:
            print e
            print("Issue opening {}\nExiting...".format(fasta))
