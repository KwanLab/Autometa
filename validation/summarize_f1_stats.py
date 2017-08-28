#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

input_path = sys.argv[1]

def summarize_f1_stats(f1_table_path):
    f1_table = pd.read_csv(f1_table_path,sep='\t')
    f1_dict = {}
    for count,row in f1_table.iterrows():
        if row['ref_genome'] != "misassembled" and str(row['ref_genome']) != "nan":
            f1_dict[row['ref_genome']] = row['F1']
    return f1_table_path,round(sum(f1_dict.values()),1),round(np.median(f1_dict.values()),1)

print summarize_f1_stats(input_path)
