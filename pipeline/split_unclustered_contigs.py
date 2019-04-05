#!/usr/bin/env python2.7

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

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
