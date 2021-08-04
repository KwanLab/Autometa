#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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

Autometa large-data-mode binning by selection of taxon sets using provided upper bound
and determined lower bound
"""


import datetime
import gzip
import logging
import os
import re
from typing import Dict, Union

import pandas as pd

from autometa.common.markers import load as load_markers
from autometa.common import kmers

from autometa.common.exceptions import TableFormatError, BinningError
from autometa.taxonomy.ncbi import NCBI
from autometa.binning.recursive_dbscan import get_clusters
from autometa.binning.utilities import (
    write_results,
    read_annotations,
    filter_taxonomy,
    reindex_bin_names,
    zero_pad_bin_names,
)


logger = logging.getLogger(__name__)


def get_kmer_embedding(
    counts: pd.DataFrame,
    cache_fpath: str,
    norm_method: str,
    pca_dimensions: int,
    embed_dimensions: int,
    embed_method: str,
) -> pd.DataFrame:
    """Retrieve kmer embeddings for provided counts by first performing kmer normalization with `norm_method`
    then PCA down to `pca_dimensions` until the normalized kmer frequencies are embedded to `embed_dimensions` using `embed_method`.

    Parameters
    ----------
    counts : pd.DataFrame
        Kmer counts where index column is 'contig' and each column is a kmer count.

    cache_fpath : str
        Path to cache embedded kmers table for later look-up/inspection.

    norm_method : str
        normalization transformation to use on kmer counts. Choices include 'am_clr', 'ilr' and 'clr'. See :func:kmers.normalize for more details.

    pca_dimensions : int
        Number of dimensions by which to initially reduce normalized kmer frequencies (Must be greater than `embed_dimensions`).

    embed_dimensions : int
        Embedding dimensions by which to reduce normalized PCA-transformed kmer frequencies (Must be less than `pca_dimensions`).

    embed_method : str
        Embedding method to use on normalized, PCA-transformed kmer frequencies. Choices include 'bhsne', 'sksne' and 'umap'. See :func:kmers.embed for more details.

    Returns
    -------
    pd.DataFrame
        [description]
    """
    # No cache dir provided so we perform normalization and embedding then return
    if not cache_fpath:
        return kmers.embed(
            kmers=kmers.normalize(counts, method=norm_method),
            embed_dimensions=embed_dimensions,
            pca_dimensions=pca_dimensions,
            method=embed_method,
        )
    # Cache was provided so we are going to first try to retrieve the cached embedding
    if os.path.exists(cache_fpath) and os.path.getsize(cache_fpath):
        # Retrieve embedding if it has already been cached and we are trying to resume.
        logger.debug(f"Found cached embeddings {cache_fpath}. Reading...")
        return pd.read_csv(cache_fpath, sep="\t", index_col="contig")
    # Cache does not exist, so we perform embedding on rank then cache
    rank_embedding = kmers.embed(
        kmers=kmers.normalize(counts, method=norm_method),
        embed_dimensions=embed_dimensions,
        pca_dimensions=pca_dimensions,
        method=embed_method,
    )
    rank_embedding.to_csv(cache_fpath, sep="\t", index=True, header=True)
    logger.debug(f"Cached embeddings to {cache_fpath}")
    return rank_embedding


def get_checkpoint_info(checkpoints_fpath: str) -> Dict[str, Union[pd.DataFrame, str]]:
    """Retrieve checkpoint information from generated binning_checkpoints.tsv

    Parameters
    ----------
    checkpoints_fpath : str
        Generated binning_checkpoints.tsv within cache directory

    Returns
    -------
    Dict[str, str, str]
        binning_checkpoints, starting canonical rank, starting rank name within starting canonical rank
        keys="binning_checkpoints", "starting_rank", "starting_rank_name_txt"
        values=pd.DataFrame, str, str
    """
    df = pd.read_csv(checkpoints_fpath, sep="\t", index_col="contig", comment="#")
    # Retrieve last column in df for most recent binning
    rank_pattern = re.compile(r"#\srank:\s(\S+)")
    rank_name_txt_pattern = re.compile(r"#\sname:\s(\S+)")
    starting_rank = None
    starting_rank_name_txt = None
    fh = (
        gzip.open(checkpoints_fpath, "rt")
        if checkpoints_fpath.endswith(".gz")
        else open(checkpoints_fpath)
    )
    for line in fh:
        rank_match = rank_pattern.search(line)
        if rank_match:
            starting_rank = rank_match.group(1)
        rank_name_txt_match = rank_name_txt_pattern.search(line)
        if rank_name_txt_match:
            starting_rank_name_txt = rank_name_txt_match.group(1)
        if starting_rank and starting_rank_name_txt:
            # At this point we have all of our variables we want to look-up
            break
    fh.close()
    logger.debug(
        f"{df.shape[1]:,} binning checkpoints found. starting at {starting_rank} with {starting_rank_name_txt}"
    )
    return {
        "binning_checkpoints": df,
        "starting_rank": starting_rank,
        "starting_rank_name_txt": starting_rank_name_txt,
    }


def checkpoint(
    checkpoints_df: pd.DataFrame,
    clustered: pd.DataFrame,
    rank: str,
    rank_name_txt: str,
    completeness: float,
    purity: float,
    domain: str,
    coverage_stddev: float,
    gc_content_stddev: float,
    cluster_method: str,
    norm_method: str,
    pca_dimensions: int,
    embed_dimensions: int,
    embed_method: str,
    min_contigs: int,
    max_partition_size: int,
    binning_checkpoints_fpath: str,
) -> pd.DataFrame:
    checkpoints_df = pd.merge(
        checkpoints_df,
        clustered[["cluster"]],
        how="outer",
        left_index=True,
        right_index=True,
    ).rename(columns={"cluster": rank_name_txt})
    checkpoints_str = checkpoints_df.to_csv(sep="\t", index=True, header=True)
    now = datetime.datetime.now()
    header = "\n".join(
        [
            f"#-- Parameters --#",
            f"# completeness: {completeness}",
            f"# purity: {purity}",
            f"# domain: {domain}",
            f"# coverage_stddev: {coverage_stddev}",
            f"# gc_content_stddev: {gc_content_stddev}",
            f"# cluster_method: {cluster_method}",
            f"# norm_method: {norm_method}",
            f"# pca_dimensions: {pca_dimensions}",
            f"# embed_dimensions: {embed_dimensions}",
            f"# embed_method: {embed_method}",
            f"# min-partition-size: {min_contigs} (max([pca_dimensions + 1, embed_dimensions + 1])",
            f"# max-partition-size: {max_partition_size}",
            f"#-- Runtime Variables --#",
            f"# rank: {rank}",
            f"# name: {rank_name_txt}",
            f"# checkpoint-shape: {checkpoints_df.shape}",
            f"# current time: {now}",
            f"# timestamp: {now.timestamp()}",
        ]
    )
    binning_checkpoints_outlines = "\n".join([header, checkpoints_str])
    if binning_checkpoints_fpath.endswith(".gz"):
        fh = gzip.open(binning_checkpoints_fpath, "wb")
        fh.write(binning_checkpoints_outlines.encode())
    else:
        fh = open(binning_checkpoints_fpath, "w")
        fh.write(binning_checkpoints_outlines)
    fh.close()
    logger.debug(
        f"Checkpoint => {rank} : {rank_name_txt} ({checkpoints_df.shape[1]:,} total checkpoints)"
    )
    return checkpoints_df


def cluster_by_taxon_partitioning(
    main: pd.DataFrame,
    counts: pd.DataFrame,
    markers: pd.DataFrame,
    norm_method: str = "am_clr",
    pca_dimensions: int = 50,
    embed_dimensions: int = 2,
    embed_method: str = "umap",
    max_partition_size: int = 10000,
    domain: str = "bacteria",
    completeness: float = 20.0,
    purity: float = 95.0,
    coverage_stddev: float = 25.0,
    gc_content_stddev: float = 5.0,
    starting_rank: str = "superkingdom",
    method: str = "dbscan",
    reverse_ranks: bool = False,
    cache: str = None,
    binning_checkpoints_fpath: str = None,
    verbose: bool = False,
) -> pd.DataFrame:
    """Perform clustering of contigs by provided `method` and use metrics to
    filter clusters that should be retained via `completeness` and `purity`
    thresholds.

    Parameters
    ----------
    main : pd.DataFrame
        index=contig,
        cols=['coverage', 'gc_content']
        taxa cols should be present if `taxonomy` is True.
        i.e. [taxid,superkingdom,phylum,class,order,family,genus,species]

    counts : pd.DataFrame
        contig kmer counts -> index_col='contig', cols=['AAAAA', 'AAAAT', ...]
        NOTE: columns will correspond to the selected k-mer count size. e.g. 3-mers would be ['AAA','AAT', ...]

    markers : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]

    domain : str, optional
        Kingdom to determine metrics (the default is 'bacteria').
        choices=['bacteria','archaea']

    completeness : float, optional
        Description of parameter `completeness` (the default is 20.).

    purity : float, optional
        purity threshold to retain cluster (the default is 95.0).
        e.g. cluster purity >= purity_cutoff

    coverage_stddev : float, optional
        cluster coverage threshold to retain cluster (the default is 25.0).

    gc_content_stddev : float, optional
        cluster GC content threshold to retain cluster (the default is 5.0).

    taxonomy : bool, optional
        Split canonical ranks and subset based on rank then attempt to find clusters (the default is True).
        taxonomic_levels = [superkingdom,phylum,class,order,family,genus,species]

    starting_rank : str, optional
        Starting canonical rank at which to begin subsetting taxonomy (the default is superkingdom).
        Choices are superkingdom, phylum, class, order, family, genus, species.

    method : str, optional
        Clustering `method` (the default is 'dbscan').
        choices = ['dbscan','hdbscan']

    reverse_ranks : bool, optional
        False - [superkingdom,phylum,class,order,family,genus,species] (Default)
        True - [species,genus,family,order,class,phylum,superkingdom]

    cache : str, optional
        Directory to cache intermediate results

    binning_checkpoints_fpath : str, optional
        File path to binning checkpoints (checkpoints are only created if the `cache` argument is provided)

    verbose : bool, optional
        log stats for each recursive_dbscan clustering iteration

    Returns
    -------
    pd.DataFrame
        main with ['cluster','completeness','purity'] columns added

    Raises
    -------
    TableFormatError
        No marker information is availble for contigs to be binned.
    FileNotFoundError
        Provided `binning_checkpoints_fpath` does not exist
    """
    if reverse_ranks:
        # species, genus, family, order, class, phylum, superkingdom
        ranks = [rank for rank in NCBI.CANONICAL_RANKS if rank != "root"]
    else:
        # superkingdom, phylum, class, order, family, genus, species
        ranks = [rank for rank in reversed(NCBI.CANONICAL_RANKS) if rank != "root"]
    # if stage is cached then we can first look to the cache before we begin subsetting main...
    clustered_contigs = set()
    num_clusters = 0
    clusters = []
    rank_embeddings = {}
    starting_rank_name_txt = None
    # Retrieve appropriate starting canonical rank and rank_name_txt from cached binning checkpoints if cache was provided
    if cache:
        if binning_checkpoints_fpath and not os.path.exists(binning_checkpoints_fpath):
            raise FileNotFoundError(binning_checkpoints_fpath)
        if not binning_checkpoints_fpath:
            binning_checkpoints_fpath = os.path.join(
                cache, "binning_checkpoints.tsv.gz"
            )
        if os.path.exists(binning_checkpoints_fpath) and os.path.getsize(
            binning_checkpoints_fpath
        ):
            checkpoint_info = get_checkpoint_info(binning_checkpoints_fpath)
            binning_checkpoints = checkpoint_info["binning_checkpoints"]
            starting_rank = checkpoint_info["starting_rank"]
            starting_rank_name_txt = checkpoint_info["starting_rank_name_txt"]
            # Update datastructures to begin at checkpoint stage.
            ## Forward fill binning annotations to most recent checkpoint and drop any contigs without bin annotations
            most_recent_binning_checkpoint = (
                binning_checkpoints.fillna(axis=1, method="ffill").iloc[:, -1].dropna()
            )
            clustered_contigs = set(
                most_recent_binning_checkpoint.index.unique().tolist()
            )
            most_recent_clustered_df = most_recent_binning_checkpoint.to_frame().rename(
                columns={starting_rank_name_txt: "cluster"}
            )
            num_clusters = most_recent_clustered_df.cluster.nunique()
            clusters.append(most_recent_clustered_df)
        else:
            logger.debug(
                f"Binning checkpoints not found. Writing checkpoints to {binning_checkpoints_fpath}"
            )
            binning_checkpoints = pd.DataFrame()

    # Subset ranks by provided (or checkpointed) starting rank
    starting_rank_index = ranks.index(starting_rank)
    ranks = ranks[starting_rank_index:]
    logger.debug(f"Using ranks: {', '.join(ranks)}")
    logger.debug(f"Max partition size set to: {max_partition_size}")
    starting_rank_name_txt_found = False
    for rank in ranks:
        # TODO: We should account for novel taxa here instead of removing 'unclassified'
        unclassified_filter = main[rank] != "unclassified"
        # group 'classified' contigs by rank
        n_taxa_in_rank = main.loc[unclassified_filter, rank].nunique()
        n_classified_contigs_in_rank = main.loc[unclassified_filter, rank].shape[0]
        logger.info(
            f"Examining {rank}: {n_taxa_in_rank:,} unique taxa ({n_classified_contigs_in_rank:,} contigs)"
        )
        # Subset counts by filtering out any unclassified according to current canonical rank and that exist in main df
        rank_counts = counts.loc[counts.index.isin(main.loc[unclassified_filter].index)]
        # cache dir structured with sub-directories corresponding to each canonical rank.
        embedding_cache_fpath = (
            os.path.join(
                cache,
                rank,
                f"{rank}.{norm_method}_pca{pca_dimensions}_{embed_method}{embed_dimensions}.tsv.gz",
            )
            if cache
            else None
        )
        # Create cache rank outdir if it does not exist
        if embedding_cache_fpath and not os.path.isdir(os.path.join(cache, rank)):
            os.makedirs(os.path.join(cache, rank))

        # Store canonical rank embedding for later lookup at lower canonical ranks
        rank_embedding = get_kmer_embedding(
            rank_counts,
            cache_fpath=embedding_cache_fpath,
            norm_method=norm_method,
            pca_dimensions=pca_dimensions,
            embed_dimensions=embed_dimensions,
            embed_method=embed_method,
        )
        rank_embeddings[rank] = rank_embedding
        # Now group by canonical rank and try embeddings specific to each name in this canonical rank
        main_grouped_by_rank = main.loc[unclassified_filter].groupby(rank)
        # Find best clusters within each rank name for the canonical rank subset
        for rank_name_txt, dff in main_grouped_by_rank:
            # Skip contig set if we are still searching for our starting taxon
            if rank_name_txt == starting_rank_name_txt:
                starting_rank_name_txt_found = True
            if starting_rank_name_txt and not starting_rank_name_txt_found:
                continue
            # Skip contig set if no contigs are classified with this rank name for this rank.
            if dff.empty:
                continue
            # Only cluster contigs that have not already been assigned a bin. (i.e. 'cluster' column value is pd.NA)
            rank_df = (
                dff.loc[dff["cluster"].isna()] if "cluster" in dff.columns else dff
            )
            # Second cluster filter with set membership to global set of clustered contigs
            if clustered_contigs:
                rank_df = rank_df.loc[~rank_df.index.isin(clustered_contigs)]
            # Post-filtering are there multiple contigs to cluster?
            if rank_df.empty:
                continue
            # First subset counts by rank_name_txt
            rank_name_txt_counts = counts.loc[counts.index.isin(rank_df.index)]
            embedding_cache_fpath = (
                os.path.join(
                    cache,
                    rank,
                    f"{rank_name_txt}.{norm_method}_pca{pca_dimensions}_{embed_method}{embed_dimensions}.tsv.gz",
                )
                if cache
                else None
            )
            # Check num. contigs for kmer embedding retrieval
            min_contigs = max([pca_dimensions + 1, embed_dimensions + 1])
            n_contigs_in_rank = rank_df.shape[0]
            # Case 1: num. contigs in rank_df is greater than the max partition size and needs to be further subset.
            if n_contigs_in_rank > max_partition_size:
                logger.debug(
                    f"{rank_name_txt} > max_partition_size ({n_contigs_in_rank:,}>{max_partition_size:,}). skipping [and caching embedding]"
                )
                rank_name_embedding = get_kmer_embedding(
                    rank_name_txt_counts,
                    cache_fpath=embedding_cache_fpath,
                    norm_method=norm_method,
                    pca_dimensions=pca_dimensions,
                    embed_dimensions=embed_dimensions,
                    embed_method=embed_method,
                )
                rank_embeddings[rank_name_txt] = rank_name_embedding
                binning_checkpoints[rank_name_txt] = pd.NA
                continue
            elif min_contigs <= n_contigs_in_rank <= max_partition_size:
                rank_name_embedding = get_kmer_embedding(
                    rank_name_txt_counts,
                    cache_fpath=embedding_cache_fpath,
                    norm_method=norm_method,
                    pca_dimensions=pca_dimensions,
                    embed_dimensions=embed_dimensions,
                    embed_method=embed_method,
                )
                logger.debug(f"{rank_name_txt} shape={rank_name_embedding.shape}")
            else:
                # Action 3: keep coordinates from higher canonical rank. i.e. kingdom embedding -> phylum embedding -> etc.
                # prev_canonical_rank_taxon_name = prev_rank_names[rank_name_txt]
                # rank_embeddings[rank_name_txt]
                rank_name_embedding = rank_embedding.loc[
                    rank_embedding.index.isin(rank_df.index)
                ]
                logger.debug(
                    f"using {rank} embedding for {rank_name_txt}. shape={rank_name_embedding.shape}"
                )
            # Add kmer embeddings into rank dataframe to use for clustering
            rank_df = pd.merge(
                rank_df,
                rank_name_embedding,
                how="left",
                left_index=True,
                right_index=True,
            )
            # Find best clusters
            logger.debug(
                f"Examining taxonomy: {rank} : {rank_name_txt} : {rank_df.shape}"
            )
            rank_name_txt_binning_df = get_clusters(
                main=rank_df,
                markers_df=markers,
                domain=domain,
                completeness=completeness,
                purity=purity,
                coverage_stddev=coverage_stddev,
                gc_content_stddev=gc_content_stddev,
                method=method,
                verbose=verbose,
            )
            # Store clustered contigs
            is_clustered = rank_name_txt_binning_df["cluster"].notnull()
            clustered = rank_name_txt_binning_df.loc[is_clustered]
            if clustered.empty:
                continue
            clustered_contigs.update({contig for contig in clustered.index})
            clustered = reindex_bin_names(df=clustered, initial_index=num_clusters)
            num_clusters += clustered.cluster.nunique()
            clusters.append(clustered)
            # Cache binning at rank_name_txt stage (rank-name-txt checkpointing)
            if cache:
                binning_checkpoints = checkpoint(
                    checkpoints_df=binning_checkpoints,
                    clustered=clustered,
                    rank=rank,
                    rank_name_txt=rank_name_txt,
                    completeness=completeness,
                    purity=purity,
                    domain=domain,
                    coverage_stddev=coverage_stddev,
                    gc_content_stddev=gc_content_stddev,
                    cluster_method=method,
                    norm_method=norm_method,
                    pca_dimensions=pca_dimensions,
                    embed_dimensions=embed_dimensions,
                    embed_method=embed_method,
                    min_contigs=min_contigs,
                    max_partition_size=max_partition_size,
                    binning_checkpoints_fpath=binning_checkpoints_fpath,
                )

    if not clusters:
        raise BinningError("Failed to recover any clusters from dataset")
    # At this point we've went through all ranks and have clusters for each canonical-rank
    # We place these into one table and then...
    clustered_df = pd.concat(clusters, sort=True)
    # ...create a dataframe for any contigs *not* in the clustered dataframe
    unclustered_df = main.loc[~main.index.isin(clustered_df.index)]
    unclustered_df["cluster"] = pd.NA
    df = pd.concat([clustered_df, unclustered_df], sort=True)
    df = zero_pad_bin_names(df)
    return df


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Autometa Large-data-mode binning by contig set selection using max-partition-size",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--kmers",
        help="Path to k-mer counts table",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--coverages",
        help="Path to metagenome coverages table",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--gc-content",
        help="Path to metagenome GC contents table",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--markers",
        help="Path to Autometa annotated markers table",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--taxonomy",
        metavar="filepath",
        help="Path to Autometa assigned taxonomies table",
        required=True,
    )
    parser.add_argument(
        "--output-binning",
        help="Path to write Autometa binning results",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--output-main",
        help="Path to write Autometa main table used during/after binning",
        metavar="filepath",
    )
    parser.add_argument(
        "--clustering-method",
        help="Clustering algorithm to use for recursive binning.",
        choices=["dbscan", "hdbscan"],
        default="dbscan",
    )
    parser.add_argument(
        "--completeness",
        help="completeness cutoff to retain cluster."
        " e.g. cluster completeness >= `completeness`",
        default=20.0,
        metavar="0 < float <= 100",
        type=float,
    )
    parser.add_argument(
        "--purity",
        help="purity cutoff to retain cluster. e.g. cluster purity >= `purity`",
        default=95.0,
        metavar="0 < float <= 100",
        type=float,
    )
    parser.add_argument(
        "--cov-stddev-limit",
        help="coverage standard deviation limit to retain cluster"
        " e.g. cluster coverage standard deviation <= `cov-stddev-limit`",
        default=25.0,
        metavar="float",
        type=float,
    )
    parser.add_argument(
        "--gc-stddev-limit",
        help="GC content standard deviation limit to retain cluster"
        " e.g. cluster GC content standard deviation <= `gc-content-stddev-limit`",
        default=5.0,
        metavar="float",
        type=float,
    )
    parser.add_argument(
        "--norm-method",
        help="kmer normalization method to use on kmer counts",
        default="am_clr",
        choices=[
            "am_clr",
            "ilr",
            "clr",
        ],
    )
    parser.add_argument(
        "--pca-dims",
        help="PCA dimensions to reduce normalized kmer frequencies prior to embedding",
        default=50,
        metavar="int",
        type=int,
    )
    parser.add_argument(
        "--embed-method",
        help="kmer embedding method to use on normalized kmer frequencies",
        default="bhsne",
        choices=[
            "bhsne",
            "umap",
            "sksne",
        ],
    )
    parser.add_argument(
        "--embed-dims",
        help="Embedding dimensions to reduce normalized kmers table after PCA.",
        default=2,
        metavar="int",
        type=int,
    )
    parser.add_argument(
        "--max-partition-size",
        help="Maximum number of contigs to consider for a recursive binning batch.",
        default=10000,
        metavar="int",
        type=int,
    )
    parser.add_argument(
        "--starting-rank",
        help="Canonical rank at which to begin subsetting taxonomy",
        default="superkingdom",
        choices=[
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
    )
    parser.add_argument(
        "--reverse-ranks",
        action="store_true",
        default=False,
        help="Reverse order at which to split taxonomy by canonical-rank."
        " When `--reverse-ranks` is given, contigs will be split in order of"
        " species, genus, family, order, class, phylum, superkingdom.",
    )
    parser.add_argument(
        "--cache",
        help="Directory to store itermediate checkpoint files during binning"
        " (If this is provided and the job fails, the script will attempt to"
        " begin from the checkpoints in this cache directory).",
        metavar="dirpath",
    )
    parser.add_argument(
        "--binning-checkpoints",
        help="File path to store itermediate contig binning results"
        " (The `--cache` argument is required for this feature). If  "
        "`--cache` is provided without this argument, a binning checkpoints file will be created.",
        metavar="filepath",
    )
    parser.add_argument(
        "--rank-filter",
        help="Taxonomy column canonical rank to subset by provided value of `--rank-name-filter`",
        default="superkingdom",
        choices=[rank for rank in NCBI.CANONICAL_RANKS if rank != "root"],
    )
    parser.add_argument(
        "--rank-name-filter",
        help="Only retrieve contigs with this name corresponding to `--rank-filter` column",
        default="bacteria",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="log debug information",
    )
    args = parser.parse_args()

    counts_df = pd.read_csv(args.kmers, sep="\t", index_col="contig")
    # First check if we are performing binning with taxonomic partitioning
    if args.taxonomy:
        main_df = read_annotations([args.coverages, args.gc_content, args.taxonomy])
        main_df = filter_taxonomy(
            df=main_df, rank=args.rank_filter, name=args.rank_name_filter
        )
    else:
        main_df = read_annotations([args.coverages, args.gc_content])
        embed_df = get_kmer_embedding(
            counts=counts_df,
            norm_method=args.norm_method,
            pca_dimensions=args.pca_dims,
            embed_dimensions=args.embed_dims,
            embed_method=args.embed_method,
            cache_fpath=None,
        )
        main_df = pd.merge(
            main_df, embed_df, how="left", left_index=True, right_index=True
        )

    # Prepare our markers dataframe
    markers_df = load_markers(args.markers, format="wide")

    # Ensure we have marker-containing contigs available to check binning quality...
    if main_df.loc[main_df.index.isin(markers_df.index)].empty:
        raise TableFormatError(
            "No markers for contigs in table. Unable to assess binning quality"
        )
    if main_df.shape[0] <= 1:
        raise BinningError("Not enough contigs in table for binning")

    logger.info(f"Selected clustering method: {args.clustering_method}")

    main_out = cluster_by_taxon_partitioning(
        main=main_df,
        counts=counts_df,
        markers=markers_df,
        norm_method=args.norm_method,
        pca_dimensions=args.pca_dims,
        embed_dimensions=args.embed_dims,
        embed_method=args.embed_method,
        max_partition_size=args.max_partition_size,
        domain=args.rank_name_filter,
        completeness=args.completeness,
        purity=args.purity,
        coverage_stddev=args.cov_stddev_limit,
        gc_content_stddev=args.gc_stddev_limit,
        starting_rank=args.starting_rank,
        method=args.clustering_method,
        reverse_ranks=args.reverse_ranks,
        cache=args.cache,
        binning_checkpoints_fpath=args.binning_checkpoints,
        verbose=args.verbose,
    )

    write_results(
        results=main_out,
        binning_output=args.output_binning,
        full_output=args.output_main,
    )


if __name__ == "__main__":
    main()
