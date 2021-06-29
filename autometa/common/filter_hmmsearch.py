#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd

from glob import glob
import os
from autometa.common.external import prodigal


def filter_markers(infpath, outfpath, cutoffs, orfs):
    for fp in [infpath, cutoffs]:
        if not os.path.exists(fp):
            raise FileNotFoundError(fp)
    if os.path.exists(outfpath) and os.path.getsize(outfpath):
        raise FileExistsError(f"{outfpath} already exists")
    col_indices = [0, 3, 4, 7]
    df = pd.read_csv(infpath, sep=r'\s+', usecols=col_indices, header=None, comment="#")
    df.rename(columns={0:"orf", 3:"sname",4:"sacc",7:"score"}, inplace=True)
    df["cleaned_sacc"] =  df["sacc"].str.replace(r'\..*$', '',regex=True)
    dff = pd.read_csv(cutoffs, sep="\t", index_col="accession")
    mdf = pd.merge(df, dff, how="left", left_on="cleaned_sacc", right_on="accession")
    mdf = mdf[mdf["score"] >= mdf["cutoff"]]
    cols = ["orf", "sacc", "sname", "score", "cutoff"]
    mdf = mdf[cols]
    if mdf.empty:
        cols = ["orf", "sacc", "sname", "score", "cutoff", "contig"]
        mdf = pd.DataFrame(columns = cols)
        mdf.set_index("contig", inplace=True)
    else:
        translations = prodigal.contigs_from_headers(orfs)    
        def translater(x):
            return translations.get(x, x.rsplit("_", 1)[0])
        mdf["contig"] = mdf["orf"].map(translater)
        mdf.set_index("contig", inplace=True)
    mdf.to_csv(outfpath, sep="\t", index=True, header=True)
    return outfpath

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Retrieves markers with provided input assembly"
    )
    parser.add_argument("--infpath", help="</path/to/domtblout>")
    parser.add_argument("--cutoffs", help="</path/to/hmm/cutoffs.tsv>")
    parser.add_argument("--markersout", help="</path/to/markers.tsv>")
    parser.add_argument("--orfs", help="</path/to/assembly.orfs.faa>")
    args = parser.parse_args()

    result = filter_markers(
        infpath=args.infpath,
        outfpath=args.markersout,
        cutoffs=args.cutoffs,
        orfs=args.orfs,
    )

if __name__ == "__main__":
    main()
