#!/usr/bin/env python

import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Split list of unclustered contigs\
    for CHTC parallelization.")
parser.add_argument('-t','--contig_tab', help='Master contig table', required=True)
parser.add_argument('-c','--cluster_column', help='Name of column for cluster', \
    default='cluster')
parser.add_argument('-b','--block_size', help='Number of contigs to place in \
    each split list file.', default=100)
parser.add_argument('-u','--unclustered_name', help='Name of unclustered group \
    in cluster column', default="unclustered")
args = vars(parser.parse_args())

###https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in xrange(0, len(l), n):
        yield l[i:i + n]
###

contig_table = pd.read_csv(args['contig_tab'], sep = "\t")
cluster_column_name = args['cluster_column']
unclustered_name = args['unclustered_name']
block_size = 10

unclustered_contig_table = contig_table.loc[contig_table[cluster_column_name] == unclustered_name]
unclustered_contig_list = unclustered_contig_table['contig'].tolist()

sublist_generator = chunks(unclustered_contig_list,block_size)
for count,sublist in enumerate(sublist_generator):
    output_name = str(count) + ".list"
    with open(output_name,"w") as outlist:
        for contig in sublist:
            outlist.write(contig + "\n")
