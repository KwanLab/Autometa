#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import logging

import pandas as pd

from autometa.common.external import prodigal
from autometa.taxonomy import majority_vote
from autometa.taxonomy.lca import LCA
from autometa.taxonomy.ncbi import NCBI_DIR
from autometa.taxonomy.ncbi import NCBI
import tempfile
import shutil
from autometa.common.metabin import MetaBin

logger = logging.getLogger(__name__)


def assign(
    method,
    outfpath,
    fasta=None,
    prot_orfs=None,
    nucl_orfs=None,
    blast=None,
    hits=None,
    lca_fpath=None,
    ncbi_dir=NCBI_DIR,
    tmpdir=None,
    usepickle=True,
    force=False,
    verbose=False,
    parallel=True,
    cpus=0,
):
    if os.path.exists(outfpath) and os.path.getsize(outfpath) and not force:
        logger.debug(
            f"FileExistsError: {outfpath}. Use force to overwrite. skipping..."
        )
        return pd.read_csv(outfpath, sep="\t", index_col="contig")
    method = method.lower()
    if method != "majority_vote":
        raise NotImplementedError(method)
    try:
        outdir = os.path.dirname(os.path.realpath(outfpath))
        tmpdir = (
            tmpdir
            if tmpdir
            else tempfile.mkdtemp(suffix=None, prefix="taxon-assignment", dir=outdir)
        )
        lca_fpath = lca_fpath if lca_fpath else os.path.join(tmpdir, "lca.tsv")
        hits = hits if hits else os.path.join(tmpdir, "hits.pkl.gz")
        blast = blast if blast else os.path.join(tmpdir, "blastp.tsv")

        def vote():
            ctg_lcas = lca.parse(lca_fpath=lca_fpath, orfs_fpath=prot_orfs)
            votes = majority_vote.rank_taxids(ctg_lcas=ctg_lcas, ncbi=lca)
            return majority_vote.write_votes(results=votes, outfpath=outfpath)

        def retrieve_lcas():
            lca.blast2lca(
                fasta=prot_orfs,
                outfpath=lca_fpath,
                blast=blast,
                hits_fpath=hits,
                force=force,
            )

        def call_orfs():
            prodigal.run(
                assembly=fasta,
                nucls_out=nucl_orfs,
                prots_out=prot_orfs,
                force=force,
                cpus=cpus,
                parallel=parallel,
            )

        logger.info(f"Assigning taxonomy via {method}. This may take a while...")

        lca = LCA(
            dbdir=ncbi_dir,
            outdir=outdir,
            usepickle=usepickle,
            verbose=verbose,
            cpus=cpus,
        )
        # Setup of taxonomy assignment sequence depending on file(s) provided
        calculation_sequence = {
            "lca_exists": [vote],
            "orfs_exists": [retrieve_lcas, vote],
            "full": [call_orfs, retrieve_lcas, vote],
        }
        # Now we need to determine which point to start the calculation...
        for fp, argname in zip(
            [lca_fpath, hits, blast, prot_orfs], ["lca", "lca", "lca", "orfs"],
        ):
            step = "full"
            if os.path.exists(fp):
                step = f"{argname}_exists"
                break

        if fasta and step == "full":
            raise ValueError(f"fasta is required if no other files are specified!")

        logger.debug(f"starting taxonomy assignment sequence from {step}")
        for calculation in calculation_sequence[step]:
            logger.debug(f"running {calculation.__name__}")
            if calculation.__name__ == "parse_bed":
                return calculation()
            calculation()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def get(fpath, assembly, ncbi_dir=NCBI_DIR, kingdom=None, outdir=None):
    if not os.path.exists(fpath) or not os.path.getsize(fpath):
        raise FileNotFoundError(fpath)
    df = pd.read_csv(fpath, sep="\t", index_col="contig")
    if df[1] <= 2:
        # fpath should only contain contig and taxid columns from voting method
        ncbi = NCBI(ncbi_dir)
        dff = ncbi.get_lineage_dataframe(df["taxid"].unique().tolist())
        df = pd.merge(
            left=df, right=dff, how="left", left_on="taxid", right_index=True,
        )
        df.to_csv(fpath, sep="\t", index=True, header=True)
        # COMBAK: Add checkpointing checksum check here
        logger.debug(f"Added canonical rank names to {fpath}")
    if "superkingdom" not in df.columns:
        raise KeyError(f"superkingdom is not in taxonomy columns {df.columns}")
    kingdoms = dict(list(df.groupby("superkingdom")))
    kingdom = kingdom.lower()
    if kingdom not in kingdoms:
        recovered_kingdoms = ", ".join(kingdoms.keys())
        raise KeyError(
            f"{kingdom} not recovered in dataset. Recovered: {recovered_kingdoms}"
        )
    bins = {}
    for superkingdom, df in kingdoms.items():
        if not kingdom:
            metabin = MetaBin(
                assembly=assembly, contigs=df.index.tolist(), outdir=outdir
            )
            bins.update({kingdom: metabin})
        if kingdom == superkingdom:
            return MetaBin(assembly=assembly, contigs=df.index.tolist(), outdir=outdir)
    return bins
