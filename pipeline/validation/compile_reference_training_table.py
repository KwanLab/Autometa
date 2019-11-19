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
