#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

Module to filter the domtbl file from hmmsearch --domtblout <filepath> using provided cutoffs
"""


import os
import logging

import pandas as pd

from autometa.common.external import prodigal


logger = logging.getLogger(__name__)


def filter_domtblout(infpath:str, cutoffs:str, orfs:str, outfpath:str = None) -> pd.DataFrame:
    #TODO::memo: Add docstring 
    for fp in [infpath, cutoffs]:
        if not os.path.exists(fp):
            raise FileNotFoundError(fp)
    if os.path.exists(outfpath) and os.path.getsize(outfpath):
        raise FileExistsError(f"{outfpath} already exists")
    col_indices = [0, 3, 4, 7]
    col_names = ["orf", "sname", "sacc", "score"]
    df = pd.read_csv(infpath, sep=r'\s+', usecols=col_indices, names=col_names, header=None, comment="#")
    # Regex: \..*$ --> search for a '.' with 0 or more trailing characters until the end of the `sacc` string
    # e.g. PF03946.9 --> PF03946
    df["cleaned_sacc"] =  df["sacc"].str.replace(r'\..*$', '', regex=True)
    logger.debug(f"{df.sacc.nunique()} unique accessions for {df.orf.nunique()} orfs")
    dff = pd.read_csv(cutoffs, sep="\t", index_col="accession")
    mdf = pd.merge(df, dff, how="left", left_on="cleaned_sacc", right_on="accession")
    mdf = mdf[mdf["score"] >= mdf["cutoff"]]
    logger.debug(f"{mdf.orf.nunique()} orfs contained {mdf.shape[0]} markers ({mdf.sacc.nunique()} unique)")
    cols = ["orf", "sacc", "sname", "score", "cutoff"]
    mdf = mdf[cols]
    # Add orf to contig mapping using header descriptions from prodigal
    if mdf.empty:
        mdf["contig"] = pd.NA
    else:
        translations = prodigal.contigs_from_headers(orfs)    
        def translater(x):
            return translations.get(x, x.rsplit("_", 1)[0])
        mdf["contig"] = mdf["orf"].map(translater)
    
    mdf.set_index("contig", inplace=True)
    if outfpath:
        mdf.to_csv(outfpath, sep="\t", index=True, header=True)
        logger.debug(f"Wrote filtered markers table to: {outfpath}")
    return mdf

def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Filters domtblout generated from hmmsearch using provided cutoffs"
    )
    parser.add_argument("--domtblout", help="Path to domtblout generated from hmmsearch -domtblout <domtblout> ... <hmmfile> <seqdb>", required=True)
    parser.add_argument("--cutoffs", help="Path to cutoffs corresponding to hmmfile used with hmmsearch <hmmfile> <seqdb>", required=True)
    parser.add_argument("--seqdb", help="Path to orfs seqdb used as input to hmmsearch ... <hmmfile> <seqdb>", required=True)
    parser.add_argument("--out", help="Path to write table of markers passing provided cutoffs", required=True)
    args = parser.parse_args()

    result = filter_domtblout(
        infpath=args.infpath,
        outfpath=args.markersout,
        cutoffs=args.cutoffs,
        orfs=args.orfs,
    )

if __name__ == "__main__":
    main()
