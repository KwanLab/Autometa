#!/usr/bin/env python

import pandas as pd
import argparse

#arguement parser
parser = argparse.ArgumentParser(description="Merge all columns of two tables \
    based on common index column. Returns table with columns sorted \
    alphanumerically")
parser.add_argument('-t1','--table_1', help='Table 1', required=True)
parser.add_argument('-t2','--table_2', help='Table 2', default='bacteria')
parser.add_argument('-i','--index_column', help='Mutual column to use as index\
    ', required=True)
parser.add_argument('-H','--how', help='Default remove redundant columns (keep\
    leftmost, with "outer"). Two keep both columns use -h "left".', \
    default='outer')
parser.add_argument('-o','--out', help='Output table name, default: \
    "merged_table.txt"', default="merged_table.txt")
args = vars(parser.parse_args())

#May need to add another argument for what to do with NA values
df1 = pd.read_csv(args['table_1'], sep = "\t")
df2 = pd.read_csv(args['table_2'], sep = "\t")

#This give you a list of columns that aren't shared by the two tables
cols_to_use = df1.columns - df2.columns
dfNew = pd.merge(df2, df1[cols_to_use], left_index=True, right_index=True, \
    how=args['how'])

#Sort columns alphanumerically
dfNew = dfNew[sorted(dfNew.columns)]
#Make sure the output table had the same number of rows as the largest input
#table.
assert(len(dfNew) == max([len(dfNew), len(df2)]))

#Write out the new dataframe to a tab-separated table
dfNew.to_csv(args['out'], sep = "\t", index = False)
