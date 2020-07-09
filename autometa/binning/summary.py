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

Script to summarize/write Autometa binning results
"""


import logging
import os
import tempfile

import pandas as pd
import numpy as np

from Bio import SeqIO
from glob import glob

from autometa import config
from autometa.taxonomy.ncbi import NCBI
from autometa.taxonomy.majority_vote import rank_taxids
from autometa.common.metabin import MetaBin


logger = logging.getLogger(__name__)


def merge_annotations(mgargs):
    binning_fpaths = {
        "bacteria": mgargs.files.bacteria_binning,
        "archaea": mgargs.files.archaea_binning,
    }
    dataframes = []
    for fpath in [mgargs.files.lengths, mgargs.files.coverages]:
        if not os.path.exists(fpath):
            raise FileNotFoundError(fpath)
        df = pd.read_csv(fpath, sep="\t", index_col="contig")
        dataframes.append(df)
    if os.path.exists(mgargs.files.taxonomy):
        df = pd.read_csv(mgargs.files.taxonomy, sep="\t", index_col="contig")
        dataframes.append(df)
    annotations = {}
    for domain, fpath in binning_fpaths.items():
        if not os.path.exists(fpath) or not os.path.getsize(fpath):
            bin_df = pd.DataFrame()
        else:
            bin_df = pd.read_csv(fpath, sep="\t", index_col="contig")
            for df in dataframes:
                bin_df = pd.merge(
                    bin_df, df, how="left", left_index=True, right_index=True
                )
        annotations.update({domain: bin_df})
    return annotations


def get_and_write_cluster_records(bin_df, mgrecords, outdir):
    """Write clusters to `outdir` given clusters `df` and metagenome `records`

    Parameters
    ----------
    mgargs : argparse.Namespace
        Metagenome args parsed from autometa.config.parse_config
    outdir : str, optional
        Path to output directory to write MetaBin fastas
        (Default is os.path.join(`mgargs.parameters.outdir`,domain))

    """
    clusters = {}
    for cluster, dff in bin_df.groupby("cluster"):
        contigs = set(dff.index)
        records = [rec for rec in mgrecords if rec.id in contigs]
        outfpath = os.path.join(outdir, f"{cluster}.fna")
        SeqIO.write(records, outfpath, "fasta")
        clusters.update({cluster: records})
    return clusters


def get_metabin_taxonomies(bin_df, ncbi_dirpath):
    """Write cluster taxonomy table to `mgargs.parameters.outdir`

    Parameters
    ----------
     mgargs : argparse.Namespace
        Metagenome args parsed from autometa.config.parse_config

    """
    canonical_ranks = NCBI.CANONICAL_RANKS
    canonical_ranks.remove("root")
    is_clustered = bin_df.cluster.notnull()
    bin_df = bin_df[is_clustered]
    tmpfpath = tempfile.mktemp()
    outcols = ["cluster", "length", "taxid", *canonical_ranks]
    bin_df[outcols].to_csv(tmpfpath, sep="\t", index=False, header=False)
    taxonomies = {}
    with open(tmpfpath) as fh:
        for line in fh:
            llist = line.strip().split("\t")
            cluster = llist[0]
            length = int(llist[1])
            taxid = int(llist[2])
            ranks = llist[3:]
            for rank, canonical_rank in zip(ranks, canonical_ranks):
                if rank != "unclassified":
                    break
            if cluster not in taxonomies:
                taxonomies.update({cluster: {canonical_rank: {taxid: length}}})
            elif canonical_rank not in taxonomies[cluster]:
                taxonomies[cluster].update({canonical_rank: {taxid: length}})
            elif taxid not in taxonomies[cluster][canonical_rank]:
                taxonomies[cluster][canonical_rank].update({taxid: length})
            else:
                taxonomies[cluster][canonical_rank][taxid] += length
    os.remove(tmpfpath)
    ncbi = NCBI(dirpath=ncbi_dirpath)
    cluster_taxonomies = rank_taxids(taxonomies, ncbi)
    cluster_taxa_df = pd.Series(data=cluster_taxonomies, name="taxid").to_frame()
    lineage_df = ncbi.get_lineage_dataframe(cluster_taxa_df.taxid.tolist(), fillna=True)
    cluster_taxa_df = pd.merge(
        cluster_taxa_df, lineage_df, how="left", left_on="taxid", right_index=True
    )
    cluster_taxa_df.index.name = "cluster"
    return cluster_taxa_df


def get_metabin_stats(bin_df, cluster_records, assembly):
    stats = []
    for cluster, dff in bin_df.groupby("cluster"):
        contigs = set(dff.index)
        seqrecords = cluster_records.get(cluster)
        metabin = MetaBin(assembly=assembly, contig_ids=contigs, seqrecords=seqrecords)
        length_weighted_coverage = np.average(
            a=dff.coverage, weights=dff.length / dff.length.sum()
        )
        stats.append(
            {
                "cluster": cluster,
                "nseqs": metabin.nseqs,
                "seqs_pct": metabin.seqs_pct,
                "size (bp)": metabin.size,
                "size_pct": metabin.size_pct,
                "N90": metabin.fragmentation_metric(quality_measure=0.9),
                "N50": metabin.fragmentation_metric(quality_measure=0.5),
                "N10": metabin.fragmentation_metric(quality_measure=0.1),
                "length_weighted_gc": metabin.length_weighted_gc,
                "length_weighted_coverage": length_weighted_coverage,
                "min_coverage": dff.coverage.min(),
                "max_coverage": dff.coverage.max(),
            }
        )
    stats_df = pd.DataFrame(stats)
    stats_df.set_index("cluster", inplace=True)
    return stats_df


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Summarize Autometa results writing genome fastas and their respective taxonomies/assembly metrics for respective metagenomes",
    )
    parser.add_argument(
        "workspace", help="Path to Autometa results workspace directory", nargs="?",
    )
    args = parser.parse_args()

    configs_search_str = os.path.join(args.workspace, "**", "metagenome_*.config")
    config_fpaths = glob(configs_search_str, recursive=True)
    for config_fpath in config_fpaths:
        mgargs = config.parse_config(config_fpath)
        annotations = merge_annotations(mgargs)
        mgrecords = [rec for rec in SeqIO.parse(mgargs.files.metagenome, "fasta")]
        for domain, bin_df in annotations.items():
            if bin_df.empty:
                logger.info(f"{domain} bin_df empty, skipping...")
                continue
            domain_outdir = os.path.join(mgargs.parameters.outdir, "metabins", domain)
            if not os.path.isdir(domain_outdir):
                os.makedirs(domain_outdir)
            logger.info(f"Writing {domain} MetaBins")
            cluster_records = get_and_write_cluster_records(
                bin_df=bin_df, mgrecords=mgrecords, outdir=domain_outdir
            )
            logger.info(f"Retrieving {domain} MetaBins' stats")
            df = get_metabin_stats(
                bin_df=bin_df,
                cluster_records=cluster_records,
                assembly=mgargs.files.metagenome,
            )
            if os.path.exists(mgargs.files.taxonomy) and os.path.getsize(
                mgargs.files.taxonomy
            ):
                logger.info(f"Retrieving {domain} MetaBins' taxonomies")
                taxa_df = get_metabin_taxonomies(
                    bin_df=bin_df, ncbi_dirpath=mgargs.databases.ncbi
                )
                df = pd.merge(df, taxa_df, left_index=True, right_index=True)

            outfpath = os.path.join(
                mgargs.parameters.outdir, "metabins", f"{domain}_metabin_summary.tsv"
            )
            df.to_csv(outfpath, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()
