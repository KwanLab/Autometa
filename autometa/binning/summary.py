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
from Bio.SeqUtils import GC
from glob import glob

from autometa import config
from autometa.taxonomy.ncbi import NCBI
from autometa.taxonomy.majority_vote import rank_taxids
from autometa.common.metabin import MetaBin
from autometa.common.markers import Markers


logger = logging.getLogger(__name__)


def get_gc_content(metagenome):
    gc_content = [
        {"contig": rec.id, "GC": GC(rec.seq)}
        for rec in SeqIO.parse(metagenome, "fasta")
    ]
    df = pd.DataFrame(gc_content).set_index("contig")
    return df


def merge_annotations(mgargs):
    """Merge all required annotations for binning summaries.

    Parameters
    ----------
    mgargs : argparse.Namespace
        metagenome args parsed from config using `config.parse_config`.

    Returns
    -------
    pd.DataFrame
        index=contig, cols=['cluster','length','coverage','taxid', *canonical_ranks]

    Raises
    ------
    FileNotFoundError
        One of lengths.tsv or coverages.tsv does not exist.
    """
    binning_fpaths = {
        "bacteria": mgargs.files.bacteria_binning,
        "archaea": mgargs.files.archaea_binning,
    }
    gc_df = get_gc_content(metagenome=mgargs.files.metagenome)
    dataframes = [gc_df]
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


def write_cluster_records(bin_df, metagenome, outdir):
    """Write clusters to `outdir` given clusters `df` and metagenome `records`

    Parameters
    ----------
    bin_df : pd.DataFrame
        Autometa binning dataframe. index='contig', cols=['cluster', ...]
    metagenome : str
        Path to metagenome fasta file
    outdir : str
        Path to output directory to write MetaBin fastas

    """
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    mgrecords = [rec for rec in SeqIO.parse(metagenome, "fasta")]
    for cluster, dff in bin_df.fillna(value={"cluster": "unclustered"}).groupby(
        "cluster"
    ):
        contigs = set(dff.index)
        records = [rec for rec in mgrecords if rec.id in contigs]
        outfpath = os.path.join(outdir, f"{cluster}.fna")
        SeqIO.write(records, outfpath, "fasta")
    return


def fragmentation_metric(df, quality_measure=0.50):
    """Describes the quality of assembled genomes that are fragmented in
        contigs of different length.

        For more information see:
            http://www.metagenomics.wiki/pdf/definition/assembly/n50

        Parameters
        ----------

        quality_measure : 0 < float < 1
            Description of parameter `quality_measure` (the default is .50).
            I.e. default measure is N50, but could use .1 for N10 or .9 for N90

        Returns
        -------
        int
            Minimum contig length to cover `quality_measure` of genome (i.e. percentile contig length)

        """
    target_size = df.length.sum() * quality_measure
    lengths = 0
    for length in df.length.sort_values(ascending=False):
        lengths += length
        if lengths > target_size:
            return length


def get_metabin_stats(bin_df, markers_fpath, assembly):
    """Retrieve statistics for all clusters recovered from Autometa binning.

    Parameters
    ----------
    bin_df : pd.DataFrame
        Autometa binning table. index=contig, cols=['cluster','length', 'GC', 'coverage', ...]
    markers_fpath : str
        </path/to/{domain}.markers.tsv
    assembly : str
        </path/to/corresponding/metagenome/assembly.fna

    Returns
    -------
    pd.DataFrame
        dataframe consisting of various metabin statistics indexed by cluster.
    """
    stats = []
    markers_df = Markers.load(markers_fpath)
    for cluster, dff in bin_df.fillna(value={"cluster": "unclustered"}).groupby(
        "cluster"
    ):
        contigs = set(dff.index)
        length_weighted_coverage = np.average(
            a=dff.coverage, weights=dff.length / dff.length.sum()
        )
        length_weighted_gc = np.average(a=dff.GC, weights=dff.length / dff.length.sum())
        num_expected_markers = markers_df.shape[1]
        pfam_counts = markers_df.loc[markers_df.index.isin(dff.index)].sum()
        if pfam_counts.empty:
            total_markers = 0
            num_single_copy_markers = 0
            num_markers_present = 0
            completeness = pd.NA
            purity = pd.NA
        else:
            total_markers = pfam_counts.sum()
            num_single_copy_markers = pfam_counts[pfam_counts == 1].count()
            num_markers_present = pfam_counts[pfam_counts >= 1].count()
            completeness = num_markers_present / num_expected_markers * 100
            purity = num_single_copy_markers / num_markers_present * 100

        stats.append(
            {
                "cluster": cluster,
                "nseqs": dff.shape[0],
                "seqs_pct": dff.shape[0] / bin_df.shape[0] * 100,
                "size (bp)": dff.length.sum(),
                "size_pct": dff.length.sum() / bin_df.length.sum() * 100,
                "N90": fragmentation_metric(dff, quality_measure=0.9),
                "N50": fragmentation_metric(dff, quality_measure=0.5),
                "N10": fragmentation_metric(dff, quality_measure=0.1),
                "length_weighted_gc": length_weighted_gc,
                "min_GC": dff.GC.min(),
                "max_GC": dff.GC.max(),
                "std_GC": dff.GC.std(),
                "length_weighted_coverage": length_weighted_coverage,
                "min_coverage": dff.coverage.min(),
                "max_coverage": dff.coverage.max(),
                "std_coverage": dff.coverage.max(),
                "completeness": completeness,
                "purity": purity,
                "num_total_markers": total_markers,
                f"num_unique_markers (expected {num_expected_markers})": num_markers_present,
                "num_single_copy_markers": num_single_copy_markers,
            }
        )
    stats_df = pd.DataFrame(stats)
    stats_df.set_index("cluster", inplace=True)
    return stats_df


def get_metabin_taxonomies(bin_df, ncbi):
    """Retrieve taxonomies of all clusters recovered from Autometa binning.

    Parameters
    ----------
    bin_df : pd.DataFrame
        Autometa binning table. index=contig, cols=['cluster','length','taxid', *canonical_ranks]
    ncbi : autometa.taxonomy.ncbi.NCBI instance
        Autometa NCBI class instance

    Returns
    -------
    pd.DataFrame
        Dataframe consisting of cluster taxonomy with taxid and canonical rank.
        Indexed by cluster
    """
    canonical_ranks = ncbi.CANONICAL_RANKS
    if "root" in canonical_ranks:
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
    cluster_taxonomies = rank_taxids(taxonomies, ncbi)
    cluster_taxa_df = pd.Series(data=cluster_taxonomies, name="taxid").to_frame()
    lineage_df = ncbi.get_lineage_dataframe(cluster_taxa_df.taxid.tolist(), fillna=True)
    cluster_taxa_df = pd.merge(
        cluster_taxa_df, lineage_df, how="left", left_on="taxid", right_index=True
    )
    cluster_taxa_df.index.name = "cluster"
    return cluster_taxa_df


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Summarize Autometa results writing genome fastas and their respective"
        " taxonomies/assembly metrics for respective metagenomes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "workspace", help="Path to Autometa results workspace directory",
    )
    parser.add_argument(
        "--write",
        help="Write autometa clustered metabins",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    configs_search_str = os.path.join(args.workspace, "**", "metagenome_*.config")
    config_fpaths = glob(configs_search_str, recursive=True)
    for config_fpath in config_fpaths:
        mgargs = config.parse_config(config_fpath)
        ncbi = NCBI(dirpath=mgargs.databases.ncbi)
        annotations = merge_annotations(mgargs)
        for domain, bin_df in annotations.items():
            if bin_df.empty:
                logger.info(f"{domain} bin_df empty, skipping...")
                continue

            if args.write:
                logger.info(f"Writing {domain} MetaBins")
                domain_outdir = os.path.join(
                    mgargs.parameters.outdir, "metabins", domain
                )
                write_cluster_records(
                    bin_df=bin_df,
                    metagenome=mgargs.files.metagenome,
                    outdir=domain_outdir,
                )
            logger.info(f"Retrieving {domain} MetaBins' stats")
            markers_fpath = (
                mgargs.files.bacteria_markers
                if domain == "bacteria"
                else mgargs.files.archaea_markers
            )
            df = get_metabin_stats(
                bin_df=bin_df,
                markers_fpath=markers_fpath,
                assembly=mgargs.files.metagenome,
            )
            if os.path.exists(mgargs.files.taxonomy) and os.path.getsize(
                mgargs.files.taxonomy
            ):
                logger.info(f"Retrieving {domain} MetaBins' taxonomies")
                taxa_df = get_metabin_taxonomies(bin_df=bin_df, ncbi=ncbi)
                df = pd.merge(
                    df, taxa_df, how="left", left_index=True, right_index=True
                )
            if args.write:
                outfpath = os.path.join(
                    mgargs.parameters.outdir,
                    "metabins",
                    f"{domain}_metabin_summary.tsv",
                )
            else:
                outfpath = os.path.join(
                    mgargs.parameters.outdir, f"{domain}_metabin_summary.tsv"
                )

            df = df.convert_dtypes()
            df.to_csv(outfpath, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()
