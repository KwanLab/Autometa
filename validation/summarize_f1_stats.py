
import pandas as pd
import numpy as np
#import os
#os.path.splitext(os.path.basename(path))[0]

input_path = "./simple_dbscan/2500Mbp/F1.tab"

def summarize_f1_stats(f1_table_path):
    f1_table = pd.read_csv(f1_table_path,sep='\t')
    f1_dict = {}
    for count,row in f1_table.iterrows():
        if row['ref_genome'] != "misassembled" and str(row['ref_genome']) != "nan":
            f1_dict[row['ref_genome']] = row['F1']
    return f1_table_path,round(sum(f1_dict.values()),1),round(np.median(f1_dict.values()),1)

#in: /home/jkwan/Sequence_data/Autometa_benchmarking
path_list = !find . -name F1.tab
for path in path_list: print summarize_f1_stats(path)
