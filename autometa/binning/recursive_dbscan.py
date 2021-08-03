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

Cluster contigs recursively searching for bins with highest completeness and purity.
"""

import gzip
import logging
import os
import re
import shutil
import tempfile
from typing import Dict, List, Tuple, Union

import pandas as pd
import numpy as np

from sklearn.cluster import DBSCAN
from hdbscan import HDBSCAN

from autometa.common.markers import load as load_markers
from autometa.common import kmers

from autometa.common.exceptions import TableFormatError, BinningError
from autometa.taxonomy.ncbi import NCBI
from autometa.binning.utilities import (
    write_results,
    read_annotations,
    add_metrics,
    apply_binning_metrics_filter,
    reindex_bin_names,
    zero_pad_bin_names,
)

pd.set_option("mode.chained_assignment", None)


logger = logging.getLogger(__name__)


def filter_taxonomy(df: pd.DataFrame, rank: str, name: str) -> pd.DataFrame:
    """Clean taxon names (by broadcasting lowercase and replacing whitespace)
    then subset by all contigs under `rank` that are equal to `name`.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing columns of canonical ranks.

    rank : str
        Canonical rank on which to apply filtering.

    name : str
        Taxon in `rank` to retrieve.

    Returns
    -------
    pd.DataFrame
        DataFrame subset by `df[rank] == name`

    Raises
    ------
    KeyError
        `rank` not in taxonomy columns.

    ValueError
        Provided `name` not found in `rank` column.
    """
    # First clean the assigned taxa by broadcasting lowercase and replacing any whitespace with underscores
    for canonical_rank in NCBI.CANONICAL_RANKS:
        if canonical_rank not in df.columns:
            continue
        df[canonical_rank] = df[canonical_rank].map(
            lambda name: name.lower()
            .replace(" ", "_")
            .replace("/", "_")
            .replace("(", "_")
            .replace(")", "_")
        )
    # Now check that the provided rank is in our dataframe
    if rank not in df.columns:
        raise KeyError(f"{rank} not in taxonomy columns: {df.columns}")
    # Perform rank filter
    filtered_df = df[df[rank] == name]
    # Check that we still have some contigs left
    if filtered_df.empty:
        raise ValueError(f"Provided name: {name} not found in {rank} column")
    logger.debug(f"{rank} filtered to {name} taxonomy. shape: {filtered_df.shape}")
    return filtered_df


def run_dbscan(
    df: pd.DataFrame,
    eps: float,
    dropcols: List[str] = [
        "cluster",
        "purity",
        "completeness",
        "coverage_stddev",
        "gc_content_stddev",
    ],
) -> pd.DataFrame:
    """Run clustering on `df` at provided `eps`.

    Notes
    -----

        * documentation for sklearn `DBSCAN <https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html>`_
        * documentation for `HDBSCAN <https://hdbscan.readthedocs.io/en/latest/index.html>`_

    Parameters
    ----------
    df : pd.DataFrame
        Contigs with embedded k-mer frequencies as ['x_1','x_2',..., 'x_ndims'] columns and 'coverage' column

    eps : float
        The maximum distance between two samples for one to be considered
        as in the neighborhood of the other. This is not a maximum bound
        on the distances of points within a cluster. This is the most
        important DBSCAN parameter to choose appropriately for your data set
        and distance function. See `DBSCAN docs <https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html>`_ for more details.

    dropcols : list, optional
        Drop columns in list from `df`
        (the default is ['cluster','purity','completeness','coverage_stddev','gc_content_stddev']).

    Returns
    -------
    pd.DataFrame
        `df` with 'cluster' column added

    Raises
    -------
    BinningError
        Dataframe is missing kmer/coverage annotations

    """
    # Ignore any errors raised from trying to drop columns that do not exist in our df.
    df.drop(columns=dropcols, inplace=True, errors="ignore")
    n_samples = df.shape[0]
    if n_samples == 1:
        clusters = pd.Series([pd.NA], index=df.index, name="cluster")
        return pd.merge(df, clusters, how="left", left_index=True, right_index=True)
    if np.any(df.isnull()):
        raise TableFormatError(
            f"df is missing {df.isnull().sum().sum()} kmer/coverage annotations"
        )
    # NOTE: all of our kmer embedded columns should correspond from "x_1" to "x_{embedding_dimensions}"
    cols = [col for col in df.columns if "x_" in col or col == "coverage"]
    # Subset what will go into clusterer to only kmer and coverage information
    X = df.loc[:, cols].to_numpy()
    # Perform clustering
    clusterer = DBSCAN(eps=eps, min_samples=1, n_jobs=-1).fit(X)
    clusters = pd.Series(clusterer.labels_, index=df.index, name="cluster")
    return pd.merge(df, clusters, how="left", left_index=True, right_index=True)


def recursive_dbscan(
    table: pd.DataFrame,
    markers_df: pd.DataFrame,
    domain: str,
    completeness_cutoff: float,
    purity_cutoff: float,
    coverage_stddev_cutoff: float,
    gc_content_stddev_cutoff: float,
    verbose: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Carry out DBSCAN, starting at eps=0.3 and continuing until there is just one
    group.

    Break conditions to speed up pipeline:
    Give up if we've got up to eps 1.3 and still no complete and pure clusters.
    Often when you start at 0.3 there are zero complete and pure clusters, because
    the groups are too small. Later, some are found as the groups enlarge enough, but
    after it becomes zero again, it is a lost cause and we may as well stop. On the
    other hand, sometimes we never find any groups, so perhaps we should give up if
    by EPS 1.3 we never find any complete/pure groups.

    Parameters
    ----------
    table : pd.DataFrame
        Contigs with embedded k-mer frequencies ('x','y'), 'coverage' and 'gc_content' columns

    markers_df : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]

    domain : str
        Kingdom to determine metrics (the default is 'bacteria').
        choices=['bacteria','archaea']

    completeness_cutoff : float
        `completeness_cutoff` threshold to retain cluster (the default is 20.0).
        e.g. cluster completeness >= completeness_cutoff

    purity_cutoff : float
        `purity_cutoff` threshold to retain cluster (the default is 95.0).
        e.g. cluster purity >= purity_cutoff

    coverage_stddev_cutoff : float
        `coverage_stddev_cutoff` threshold to retain cluster (the default is 25.0).
        e.g. cluster coverage std.dev. <= coverage_stddev_cutoff

    gc_content_stddev_cutoff : float
        `gc_content_stddev_cutoff` threshold to retain cluster (the default is 5.0).
        e.g. cluster gc_content std.dev. <= gc_content_stddev_cutoff

    verbose : bool
        log stats for each recursive_dbscan clustering iteration.

    Returns
    -------
    2-tuple
        (pd.DataFrame(<passed cutoffs>), pd.DataFrame(<failed cutoffs>))
        DataFrames consisting of contigs that passed/failed clustering cutoffs, respectively.

        DataFrame:
            index = contig
            columns = ['x,'y','coverage','gc_content','cluster','purity','completeness','coverage_stddev','gc_content_stddev']
    """
    eps = 0.3
    step = 0.1
    clustering_rounds = {0: 0}
    n_clusters = float("inf")
    best_median = float("-inf")
    best_df = pd.DataFrame()
    while n_clusters > 1:
        binned_df = run_dbscan(table, eps)
        df, metrics_df = add_metrics(df=binned_df, markers_df=markers_df, domain=domain)
        filtered_df = apply_binning_metrics_filter(
            df=metrics_df,
            completeness_cutoff=completeness_cutoff,
            purity_cutoff=purity_cutoff,
            coverage_stddev_cutoff=coverage_stddev_cutoff,
            gc_content_stddev_cutoff=gc_content_stddev_cutoff,
        )
        median_completeness = filtered_df.completeness.median()
        if median_completeness >= best_median:
            best_median = median_completeness
            best_df = df
        # Count the number of clusters
        n_clusters = df["cluster"].nunique()
        if n_clusters in clustering_rounds:
            clustering_rounds[n_clusters] += 1
        else:
            clustering_rounds[n_clusters] = 1
        # We speed this up if we are getting a lot of tables with the same number of clusters
        if clustering_rounds[n_clusters] > 10:
            step *= 10
        if median_completeness == 0:
            clustering_rounds[0] += 1
        if clustering_rounds[0] >= 10:
            break
        if verbose:
            logger.debug(
                f"EPS: {eps:3.2f} Clusters: {n_clusters:,}"
                f" Completeness: Median={median_completeness:4.2f} Best={best_median:4.2f}"
            )
        eps += step
    if best_df.empty:
        if verbose:
            logger.debug("No complete or pure clusters found")
        return pd.DataFrame(), table

    clustered_df = apply_binning_metrics_filter(
        df=best_df,
        completeness_cutoff=completeness_cutoff,
        purity_cutoff=purity_cutoff,
        coverage_stddev_cutoff=coverage_stddev_cutoff,
        gc_content_stddev_cutoff=gc_content_stddev_cutoff,
    )
    unclustered_df = best_df.loc[~best_df.index.isin(clustered_df.index)]
    if verbose:
        logger.debug(f"Best completeness median: {best_median:4.2f}")
    logger.debug(
        f"clustered: {clustered_df.shape[0]:,} unclustered: {unclustered_df.shape[0]:,}"
    )
    return clustered_df, unclustered_df


def run_hdbscan(
    df: pd.DataFrame,
    min_cluster_size: int,
    min_samples: int,
    cache_dir: str = None,
    dropcols: List[str] = [
        "cluster",
        "purity",
        "completeness",
        "coverage_stddev",
        "gc_content_stddev",
    ],
) -> pd.DataFrame:
    """Run clustering on `df` at provided `min_cluster_size`.

    Notes
    -----

        * reasoning for parameter: `cluster_selection_method <https://hdbscan.readthedocs.io/en/latest/parameter_selection.html#leaf-clustering>`_
        * reasoning for parameters: `min_cluster_size and min_samples <https://hdbscan.readthedocs.io/en/latest/parameter_selection.html>`_
        * documentation for `HDBSCAN <https://hdbscan.readthedocs.io/en/latest/index.html>`_

    Parameters
    ----------
    df : pd.DataFrame
        Contigs with embedded k-mer frequencies as ['x','y'] columns and optionally 'coverage' column

    min_cluster_size : int
        The minimum size of clusters; single linkage splits that contain
        fewer points than this will be considered points "falling out" of a
        cluster rather than a cluster splitting into two new clusters.

    min_samples : int
        The number of samples in a neighborhood for a point to be
        considered a core point.

    cache_dir : str, optional
        Used to cache the output of the computation of the tree.
        By default, no caching is done. If a string is given, it is the
        path to the caching directory.

    dropcols : list, optional
        Drop columns in list from `df`
        (the default is ['cluster','purity','completeness']).

    Returns
    -------
    pd.DataFrame
        `df` with 'cluster' column added

    Raises
    -------
    ValueError
        sets `usecols` and `dropcols` may not share elements
    TableFormatError
        `df` is missing k-mer or coverage annotations.

    """
    # Ignore any errors raised from trying to drop columns that do not exist in our df.
    df.drop(columns=dropcols, inplace=True, errors="ignore")
    n_samples = df.shape[0]
    if n_samples == 1:
        clusters = pd.Series([pd.NA], index=df.index, name="cluster")
        return pd.merge(df, clusters, how="left", left_index=True, right_index=True)
    if np.any(df.isnull()):
        raise TableFormatError(
            f"df is missing {df.isnull().sum().sum()} kmer/coverage annotations"
        )
    # NOTE: all of our kmer embedded columns should correspond from "x_1" to "x_{embedding_dimensions}"
    cols = [col for col in df.columns if "x_" in col or col == "coverage"]
    # Subset what will go into clusterer to only kmer and coverage information
    X = df.loc[:, cols].to_numpy()
    # Perform clustering
    clusterer = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method="leaf",
        allow_single_cluster=True,
        memory=cache_dir,
    ).fit(X)
    clusters = pd.Series(clusterer.labels_, index=df.index, name="cluster")
    return pd.merge(df, clusters, how="left", left_index=True, right_index=True)


def recursive_hdbscan(
    table: pd.DataFrame,
    markers_df: pd.DataFrame,
    domain: str,
    completeness_cutoff: float,
    purity_cutoff: float,
    coverage_stddev_cutoff: float,
    gc_content_stddev_cutoff: float,
    verbose: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Recursively run HDBSCAN starting with defaults and iterating the min_samples
     and min_cluster_size until only 1 cluster is recovered.

    Parameters
    ----------
    table : pd.DataFrame
        Contigs with embedded k-mer frequencies ('x','y'), 'coverage' and 'gc_content' columns

    markers_df : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]

    domain : str
        Kingdom to determine metrics (the default is 'bacteria').
        choices=['bacteria','archaea']

    completeness_cutoff : float
        `completeness_cutoff` threshold to retain cluster.
        e.g. cluster completeness >= completeness_cutoff

    purity_cutoff : float
        `purity_cutoff` threshold to retain cluster.
        e.g. cluster purity >= purity_cutoff

    coverage_stddev_cutoff : float
        `coverage_stddev_cutoff` threshold to retain cluster.
        e.g. cluster coverage std.dev. <= coverage_stddev_cutoff

    gc_content_stddev_cutoff : float
        `gc_content_stddev_cutoff` threshold to retain cluster.
        e.g. cluster gc_content std.dev. <= gc_content_stddev_cutoff

    verbose : bool
        log stats for each recursive_dbscan clustering iteration.

    Returns
    -------
    2-tuple
        (pd.DataFrame(<passed cutoffs>), pd.DataFrame(<failed cutoffs>))
        DataFrames consisting of contigs that passed/failed clustering
        cutoffs, respectively.

        DataFrame:
            index = contig
            columns = ['x_1','x_2','coverage','gc_content','cluster','purity','completeness','coverage_stddev','gc_content_stddev']
    """
    max_min_cluster_size = 10000
    max_min_samples = 10
    min_cluster_size = 2
    min_samples = 1
    n_clusters = float("inf")
    best_median = float("-inf")
    best_df = pd.DataFrame()
    cache_dir = tempfile.mkdtemp()
    while n_clusters > 1:
        binned_df = run_hdbscan(
            table,
            min_cluster_size=min_cluster_size,
            min_samples=min_samples,
            cache_dir=cache_dir,
        )
        df, metrics_df = add_metrics(df=binned_df, markers_df=markers_df, domain=domain)
        filtered_df = apply_binning_metrics_filter(
            df=metrics_df,
            completeness_cutoff=completeness_cutoff,
            purity_cutoff=purity_cutoff,
            coverage_stddev_cutoff=coverage_stddev_cutoff,
            gc_content_stddev_cutoff=gc_content_stddev_cutoff,
        )
        median_completeness = filtered_df.completeness.median()
        if median_completeness >= best_median:
            best_median = median_completeness
            best_df = df

        n_clusters = df["cluster"].nunique()

        if verbose:
            logger.debug(
                f"(min_samples, min_cluster_size): ({min_samples}, {min_cluster_size}) clusters: {n_clusters}"
                f" median completeness (current, best): ({median_completeness:4.2f}, {best_median:4.2f})"
            )

        if min_cluster_size >= max_min_cluster_size:
            shutil.rmtree(cache_dir)
            cache_dir = tempfile.mkdtemp()
            min_samples += 1
            min_cluster_size = 2
        else:
            min_cluster_size += 10

        if filtered_df.empty:
            min_cluster_size += 100

        if min_samples >= max_min_samples:
            max_min_cluster_size *= 2

    # clean up cache now that we are out of while loop
    shutil.rmtree(cache_dir)
    # Check our df is not empty from while loop
    if best_df.empty:
        if verbose:
            logger.debug("No complete or pure clusters found")
        return pd.DataFrame(), table

    # Now split our clustered/unclustered contigs into two dataframes
    # First keep only clusters passing all binning filters from clustering dataframe
    clustered_df = apply_binning_metrics_filter(
        df=best_df,
        completeness_cutoff=completeness_cutoff,
        purity_cutoff=purity_cutoff,
        coverage_stddev_cutoff=coverage_stddev_cutoff,
        gc_content_stddev_cutoff=gc_content_stddev_cutoff,
    )
    # Second, using the clustered contigs locate contigs that were *not* clustered (Note the use of the ~ operator)
    unclustered_df = best_df.loc[~best_df.index.isin(clustered_df.index)]
    if verbose:
        logger.debug(f"Best completeness median: {best_median:4.2f}")
    logger.debug(
        f"clustered: {clustered_df.shape[0]} unclustered: {unclustered_df.shape[0]}"
    )
    return clustered_df, unclustered_df


def get_clusters(
    main: pd.DataFrame,
    markers_df: pd.DataFrame,
    domain: str,
    completeness: float,
    purity: float,
    coverage_stddev: float,
    gc_content_stddev: float,
    method: str,
    verbose: bool = False,
) -> pd.DataFrame:
    """Find best clusters retained after applying metrics filters.

    Parameters
    ----------
    main : pd.DataFrame
        index=contig,
        cols=['x','y','coverage','gc_content']

    markers_df : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]

    domain : str
        Kingdom to determine metrics.
        choices=['bacteria','archaea'].

    completeness : float
        completeness threshold to retain cluster.
        e.g. cluster completeness >= completeness

    purity : float
        `purity` threshold to retain cluster.
        e.g. cluster purity >= purity

    coverage_stddev : float
        cluster coverage std.dev. threshold to retain cluster.
        e.g. cluster coverage std.dev. <= coverage_stddev

    gc_content_stddev : float
        cluster GC content std.dev. threshold to retain cluster.
        e.g. cluster GC content std.dev. <= gc_content_stddev

    method : str
        Description of parameter `method`.
        choices = ['dbscan','hdbscan']

    verbose : bool
        log stats for each recursive_dbscan clustering iteration

    Returns
    -------
    pd.DataFrame
        `main` with ['cluster','completeness','purity','coverage_stddev','gc_content_stddev'] columns added

    """
    num_clusters = 0
    clusters = []
    recursive_clusterers = {"dbscan": recursive_dbscan, "hdbscan": recursive_hdbscan}
    if method not in recursive_clusterers:
        raise ValueError(
            f"Method: {method} not a choice. choose between {recursive_clusterers.keys()}"
        )
    clusterer = recursive_clusterers[method]

    # Continue until clusters are no longer being recovered
    # break when either clustered_df or unclustered_df is empty
    while True:
        clustered_df, unclustered_df = clusterer(
            main,
            markers_df,
            domain,
            completeness,
            purity,
            coverage_stddev,
            gc_content_stddev,
            verbose=verbose,
        )
        # No contigs can be clustered, label as unclustered and add the final df
        # of (unclustered) contigs
        if clustered_df.empty:
            unclustered_df = unclustered_df.assign(cluster=pd.NA)
            clusters.append(unclustered_df)
            break

        clustered_df = reindex_bin_names(df=clustered_df, initial_index=num_clusters)

        # All contigs have now been clustered, add the final df of (clustered) contigs
        if unclustered_df.empty:
            clusters.append(clustered_df)
            break
        # Store clustered contigs
        num_clusters += clustered_df.cluster.nunique()
        clusters.append(clustered_df)
        # continue with unclustered contigs
        main = unclustered_df
    # concatenate all clustering rounds in clusters list and sort by bin number
    df = pd.concat(clusters, sort=True)
    for metric in [
        "purity",
        "completeness",
        "coverage_stddev",
        "gc_content_stddev",
    ]:
        # Special case where no clusters are found in first round.
        # i.e. cluster = pd.NA for all contigs
        if metric not in df.columns:
            df[metric] = pd.NA

    df = zero_pad_bin_names(df)
    return df


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
            # Now check num. contigs for kmer embedding retrieval
            # We are accounting for three cases.
            # Case 1: num. contigs in rank_df is greater than the max partition size and needs to be further subset.
            #   Action 1: embed these contigs then continue to recurse s.t. their higher canonical rank coords.
            #   may be looked up at the lower canonical rank (See Case 2)
            # Case 2: num. contigs in rank_df is less than max partition size and greater than min contigs
            #   Action 2: perform normalization and embedding on counts for specific ranks' contigs
            # Case 3: num. contigs in rank_df is less than the minimum number of contigs required for embedding
            #   Action 3: keep coordinates from higher canonical rank. i.e. kingdom embedding -> phylum embedding -> etc.
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
                binning_checkpoints = pd.merge(
                    binning_checkpoints,
                    clustered[["cluster"]],
                    how="outer",
                    left_index=True,
                    right_index=True,
                ).rename(columns={"cluster": rank_name_txt})
                binning_checkpoints_str = binning_checkpoints.to_csv(
                    sep="\t", index=True, header=True
                )
                header = "\n".join(
                    [
                        f"#-- Parameters --#",
                        f"# completeness: {completeness}",
                        f"# purity: {purity}",
                        f"# domain: {domain}",
                        f"# coverage_stddev: {coverage_stddev}",
                        f"# gc_content_stddev: {gc_content_stddev}",
                        f"# method: {method}",
                        f"# norm_method: {norm_method}",
                        f"# pca_dimensions: {pca_dimensions}",
                        f"# embed_dimensions: {embed_dimensions}",
                        f"# embed_method: {embed_method}",
                        f"# min-partition-size: {min_contigs} (max([pca_dimensions + 1, embed_dimensions + 1])",
                        f"# max-partition-size: {max_partition_size}",
                        f"#-- Runtime Variables --#",
                        f"# rank: {rank}",
                        f"# name: {rank_name_txt}",
                        f"# checkpoint-shape: {binning_checkpoints.shape}",
                    ]
                )
                binning_checkpoints_outlines = "\n".join(
                    [header, binning_checkpoints_str]
                )

                if binning_checkpoints_fpath.endswith(".gz"):
                    fh = gzip.open(binning_checkpoints_fpath, "wb")
                    fh.write(binning_checkpoints_outlines.encode())
                else:
                    fh = open(binning_checkpoints_fpath, "w")
                    fh.write(binning_checkpoints_outlines)
                fh.close()
                logger.debug(
                    f"Checkpoint => {rank} : {rank_name_txt} ({binning_checkpoints.shape[1]:,} total checkpoints)"
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
        description="Perform marker gene guided binning of "
        "metagenome contigs using annotations (when available) of sequence "
        "composition, coverage and homology.",
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
        "--taxonomy",
        metavar="filepath",
        help="Path to Autometa assigned taxonomies table",
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

    # Perform clustering w/o taxonomy
    if not args.taxonomy:
        main_out = get_clusters(
            main=main_df,
            markers_df=markers_df,
            domain=args.rank_name_filter,
            completeness=args.completeness,
            purity=args.purity,
            coverage_stddev=args.coverage_stddev,
            gc_content_stddev=args.gc_content_stddev,
            method=args.clustering_method,
            verbose=args.verbose,
        )
    else:
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

    # reindex_bin_names(df=, cluster_col=, initial_index=)

    write_results(
        results=main_out,
        binning_output=args.output_binning,
        full_output=args.output_main,
    )


if __name__ == "__main__":
    main()
