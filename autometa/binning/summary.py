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


def write_clusters(mgargs):
    """Write clusters to `outdir` given clusters `df` and metagenome `records`

    Parameters
    ----------
    mgargs : argparse.Namespace
        Metagenome args parsed from autometa.config.parse_config
    outdir : str, optional
        Path to output directory to write MetaBin fastas
        (Default is os.path.join(`mgargs.parameters.outdir`,domain))

    """
    records = [rec for rec in SeqIO.parse(mgargs.files.metagenome, "fasta")]
    binning_fpaths = {
        "bacteria": mgargs.files.bacteria_binning,
        "archaea": mgargs.files.archaea_binning,
    }
    for domain, binning_fpath in binning_fpaths.items():
        if not os.path.exists(binning_fpath):
            logger.info(f"{binning_fpath} not found, skipping...")
            continue
        df = pd.read_csv(binning_fpath, sep="\t", index_col="contig")
        outdir = os.path.join(mgargs.parameters.outdir, "metabins", domain)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        for cluster, dff in df.groupby("cluster"):
            contigs = set(dff.index)
            recs = [rec for rec in records if rec.id in contigs]
            outfpath = os.path.join(outdir, f"{cluster}.fna")
            SeqIO.write(recs, outfpath, "fasta")


def get_metabin_taxonomies(mgargs):
    """Write cluster taxonomy table to `mgargs.parameters.outdir`

    Parameters
    ----------
     mgargs : argparse.Namespace
        Metagenome args parsed from autometa.config.parse_config

    """
    binning_fpaths = {
        "bacteria": mgargs.files.bacteria_binning,
        "archaea": mgargs.files.archaea_binning,
    }
    taxa_df = pd.read_csv(mgargs.files.taxonomy, sep="\t", index_col="contig")
    length_df = pd.read_csv(mgargs.files.lengths, sep="\t", index_col="contig")
    ncbi = NCBI(dirpath=mgargs.databases.ncbi)
    canonical_ranks = NCBI.CANONICAL_RANKS
    canonical_ranks.remove("root")
    domain_taxonomies = []
    for domain, binning_fpath in binning_fpaths.items():
        if not os.path.exists(binning_fpath):
            logger.info(f"{binning_fpath} not found, skipping...")
            continue
        bin_df = pd.read_csv(binning_fpath, sep="\t", index_col="contig")
        is_clustered = bin_df.cluster.notnull()
        bin_df = bin_df[is_clustered]
        for df in [length_df, taxa_df]:
            bin_df = pd.merge(bin_df, df, how="left", left_index=True, right_index=True)
        tmpfpath = tempfile.mktemp()
        outcols = ["cluster", "length", "taxid"]
        bin_df[outcols].to_csv(tmpfpath, sep="\t", index=False, header=False)
        taxonomies = {}
        with open(tmpfpath) as fh:
            for line in fh:
                cluster, length, taxid = line.strip().split("\t")
                rank = "root"
                for canonical_rank in canonical_ranks:
                    if canonical_rank != "unclassified":
                        rank = canonical_rank
                        break
                taxid = int(taxid)
                length = int(length)
                if cluster not in taxonomies:
                    taxonomies.update({cluster: {rank: {taxid: length}}})
                elif rank not in taxonomies[cluster]:
                    taxonomies[cluster].update({rank: {taxid: length}})
                elif taxid not in taxonomies[cluster][rank]:
                    taxonomies[cluster][rank].update({taxid: length})
                else:
                    taxonomies[cluster][rank][taxid] += length
        os.remove(tmpfpath)
        cluster_taxonomies = rank_taxids(taxonomies, ncbi)
        cluster_taxa_df = pd.Series(data=cluster_taxonomies, name="taxid").to_frame()
        lineage_df = ncbi.get_lineage_dataframe(
            cluster_taxa_df.taxid.tolist(), fillna=True
        )
        cluster_taxa_df = pd.merge(
            cluster_taxa_df, lineage_df, how="left", left_on="taxid", right_index=True
        )
        cluster_taxa_df.index.name = "cluster"
        cluster_taxa_df.reset_index(inplace=True)
        domain_taxonomies.append(cluster_taxa_df)

    master_taxa_df = pd.concat(domain_taxonomies, ignore_index=True)
    outcols = [*reversed(canonical_ranks), "taxid"]
    outcols.remove("superkingdom")
    # outfpath = os.path.join(
    #     mgargs.parameters.outdir, "metabins", "cluster_taxonomy.tsv"
    # )
    # master_taxa_df[outcols].to_csv(outfpath, sep="\t", index=False, header=True)
    master_taxa_df.set_index(["superkingdom", "cluster"], inplace=True)
    return master_taxa_df[outcols]


def get_metabin_stats(mgargs):
    cov_df = pd.read_csv(mgargs.files.coverages, sep="\t", index_col="contig")
    length_df = pd.read_csv(mgargs.files.lengths, sep="\t", index_col="contig")
    binning_fpaths = {
        "bacteria": mgargs.files.bacteria_binning,
        "archaea": mgargs.files.archaea_binning,
    }
    stats = []
    for domain, binning_fpath in binning_fpaths.items():
        if not os.path.exists(binning_fpath):
            logger.info(f"{binning_fpath} not found, skipping...")
            continue
        bin_df = pd.read_csv(binning_fpath, sep="\t", index_col="contig")
        for df in [length_df, cov_df]:
            bin_df = pd.merge(bin_df, df, how="left", left_index=True, right_index=True)
        for cluster, dff in bin_df.groupby("cluster"):
            contigs = set(dff.index)
            metabin = MetaBin(
                assembly=mgargs.files.metagenome,
                contigs=contigs,
                outdir=mgargs.parameters.outdir,
            )
            length_weighted_coverage = np.average(
                a=dff.coverage, weights=dff.length / dff.length.sum()
            )
            stats.append(
                {
                    "superkingdom": domain,
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
    stats_df.set_index(["superkingdom", "cluster"], inplace=True)
    # outfpath = os.path.join(mgargs.parameters.outdir, "metabins", "cluster_stats.tsv")
    # stats_df.to_csv(outfpath, sep="\t", index=True, header=True)
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

        # df = merge_annotations(mgargs)
        logger.info("Writing MetaBins")
        write_clusters(mgargs)
        logger.info("Retrieving MetaBins' taxonomies")
        taxa_df = get_metabin_taxonomies(mgargs)
        logger.info("Retrieving MetaBins' stats")
        stats_df = get_metabin_stats(mgargs)

        df = pd.merge(stats_df, taxa_df, left_index=True, right_index=True)
        outfpath = os.path.join(
            mgargs.parameters.outdir, "metabins", "metabin_summary.tsv"
        )
        df.to_csv(outfpath, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()
