#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Iza1.3ak Miller, Jason C. Kwan

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

import logging
import shutil
import tempfile
from typing import List, Tuple

import pandas as pd
import numpy as np

from sklearn.cluster import DBSCAN
from hdbscan import HDBSCAN

from autometa.common.markers import load as load_markers
from autometa.common import kmers

from autometa.common.exceptions import TableFormatError, BinningError
from autometa.taxonomy.ncbi import NCBI

pd.set_option("mode.chained_assignment", None)


logger = logging.getLogger(__name__)


def add_metrics(
    df: pd.DataFrame, markers_df: pd.DataFrame, domain: str = "bacteria"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Adds the completeness and purity metrics to each respective contig in df.

    Parameters
    ----------
    df : pd.DataFrame
        index='contig' cols=['x','y','coverage','gc_content','cluster']
    markers_df : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]
    domain : str, optional
        Kingdom to determine metrics (the default is 'bacteria').
        choices=['bacteria','archaea']

    Returns
    -------
    2-tuple
        `df` with added cols=['completeness', 'purity']
        pd.DataFrame(index=clusters,  cols=['completeness', 'purity'])

    Raises
    -------
    KeyError
        `domain` is not "bacteria" or "archaea"

    """
    domain = domain.lower()
    marker_sets = {"bacteria": 139.0, "archaea": 162.0}
    if domain not in marker_sets:
        raise KeyError(f"{domain} is not bacteria or archaea!")
    expected_number = marker_sets[domain]
    if "cluster" in df.columns:
        # join cluster and marker data, group by cluster
        temp = df.join(markers_df, how="outer").groupby("cluster")
        # count present
        nunique_markers = temp[list(markers_df.columns)].sum().ge(1).sum(axis=1)
        # count single copy
        num_single_copy_markers = temp[list(markers_df.columns)].sum().eq(1).sum(axis=1)
        # calculate completeness/purity
        completeness = nunique_markers / expected_number * 100
        purity = num_single_copy_markers / nunique_markers * 100
        coverage_stddev = temp["coverage"].std()
        gc_content_stddev = temp["gc_content"].std()
        completeness = completeness.to_frame()
        purity = purity.to_frame()
        coverage_stddev = coverage_stddev.to_frame()
        gc_content_stddev = gc_content_stddev.to_frame()
        metrics_df = pd.concat(
            [completeness, purity, coverage_stddev, gc_content_stddev], axis=1
        )
        metrics_df.columns = [
            "completeness",
            "purity",
            "coverage_stddev",
            "gc_content_stddev",
        ]
        merged_df = pd.merge(df, metrics_df, left_on="cluster", right_index=True)
    # Account for exceptions where clusters were not recovered
    else:
        metrics_df = pd.DataFrame(
            [
                {
                    "contig": contig,
                    "cluster": pd.NA,
                    "completeness": pd.NA,
                    "purity": pd.NA,
                    "coverage_stddev": pd.NA,
                    "gc_content_stddev": pd.NA,
                }
                for contig in df.index
            ]
        )
        metric_cols = ["completeness", "purity", "coverage_stddev", "gc_content_stddev"]
        merged_df = df.copy()
        for metric in metric_cols:
            merged_df[metric] = pd.NA
    return merged_df, metrics_df


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
        (the default is ['cluster','purity','completeness']).

    Returns
    -------
    pd.DataFrame
        `df` with 'cluster' column added

    Raises
    -------
    BinningError
    ------------
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
        DataFrames consisting of contigs that passed/failed clustering
        cutoffs, respectively.

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
        completeness_filter = metrics_df["completeness"] >= completeness_cutoff
        purity_filter = metrics_df["purity"] >= purity_cutoff
        coverage_filter = metrics_df["coverage_stddev"] <= coverage_stddev_cutoff
        gc_content_filter = metrics_df["gc_content_stddev"] <= gc_content_stddev_cutoff
        filtered_df = metrics_df[
            completeness_filter & purity_filter & coverage_filter & gc_content_filter
        ]
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

    completeness_filter = best_df["completeness"] >= completeness_cutoff
    purity_filter = best_df["purity"] >= purity_cutoff
    coverage_filter = best_df["coverage_stddev"] <= coverage_stddev_cutoff
    gc_content_filter = best_df["gc_content_stddev"] <= gc_content_stddev_cutoff
    clustered_df = best_df[
        completeness_filter & purity_filter & coverage_filter & gc_content_filter
    ]
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
        `completeness_cutoff` threshold to retain cluster (the default is 20.0).
        e.g. cluster completeness >= completeness_cutoff

    purity_cutoff : float
        `purity_cutoff` threshold to retain cluster (the default is 95.00).
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
        DataFrames consisting of contigs that passed/failed clustering
        cutoffs, respectively.

        DataFrame:
            index = contig
            columns = ['x,'y','coverage','gc_content','cluster','purity','completeness','coverage_stddev','gc_content_stddev']
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
        completeness_filter = metrics_df["completeness"] >= completeness_cutoff
        purity_filter = metrics_df["purity"] >= purity_cutoff
        coverage_filter = metrics_df["coverage_stddev"] <= coverage_stddev_cutoff
        gc_content_filter = metrics_df["gc_content_stddev"] <= gc_content_stddev_cutoff
        filtered_df = metrics_df[
            completeness_filter & purity_filter & coverage_filter & gc_content_filter
        ]
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

        if metrics_df[completeness_filter & purity_filter].empty:
            min_cluster_size += 100

        if min_samples >= max_min_samples:
            max_min_cluster_size *= 2

    shutil.rmtree(cache_dir)

    if best_df.empty:
        if verbose:
            logger.debug("No complete or pure clusters found")
        return pd.DataFrame(), table

    completeness_filter = best_df["completeness"] >= completeness_cutoff
    purity_filter = best_df["purity"] >= purity_cutoff
    coverage_filter = best_df["coverage_stddev"] <= coverage_stddev_cutoff
    gc_content_filter = best_df["gc_content_stddev"] <= gc_content_stddev_cutoff
    clustered_df = best_df[
        completeness_filter & purity_filter & coverage_filter & gc_content_filter
    ]
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
    """Find best clusters retained after applying `completeness` and `purity` filters.

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
        `main` with ['cluster','completeness','purity'] columns added
    """
    num_clusters = 0
    clusters = []
    recursive_clusterers = {"dbscan": recursive_dbscan, "hdbscan": recursive_hdbscan}
    if method not in recursive_clusterers:
        raise ValueError(f"Method: {method} not a choice. choose b/w dbscan & hdbscan")
    clusterer = recursive_clusterers[method]

    # Continue while unclustered are remaining
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

        translation = {
            c: f"bin_{1+i+num_clusters:04d}"
            for i, c in enumerate(clustered_df.cluster.unique())
        }

        def rename_cluster(c):
            return translation[c]

        clustered_df.cluster = clustered_df.cluster.map(rename_cluster)

        # All contigs have now been clustered, add the final df of (clustered) contigs
        if unclustered_df.empty:
            clusters.append(clustered_df)
            break
        # Store clustered contigs
        num_clusters += clustered_df.cluster.nunique()
        clusters.append(clustered_df)
        # continue with unclustered contigs
        main = unclustered_df
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
    return df


def binning(
    main: pd.DataFrame,
    markers: pd.DataFrame,
    domain: str = "bacteria",
    completeness: float = 20.0,
    purity: float = 95.0,
    coverage_stddev: float = 25.0,
    gc_content_stddev: float = 5.0,
    taxonomy: bool = True,
    starting_rank: str = "superkingdom",
    method: str = "dbscan",
    reverse_ranks: bool = False,
    verbose: bool = False,
) -> pd.DataFrame:
    """Perform clustering of contigs by provided `method` and use metrics to
    filter clusters that should be retained via `completeness` and `purity`
    thresholds.

    Parameters
    ----------
    main : pd.DataFrame
        index=contig,
        cols=['x','y', 'coverage', 'gc_content']
        taxa cols should be present if `taxonomy` is True.
        i.e. [taxid,superkingdom,phylum,class,order,family,genus,species]

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
    """
    # First check needs to ensure we have markers available to check binning quality...
    if main.loc[main.index.isin(markers.index)].empty:
        raise TableFormatError(
            "No markers for contigs in table. Unable to assess binning quality"
        )
    if main.shape[0] <= 1:
        raise BinningError("Not enough contigs in table for binning")

    logger.info(f"Using {method} clustering method")
    if not taxonomy:
        return get_clusters(
            main=main,
            markers_df=markers,
            domain=domain,
            completeness=completeness,
            purity=purity,
            coverage_stddev=coverage_stddev,
            gc_content_stddev=gc_content_stddev,
            method=method,
            verbose=verbose,
        )

    # Use taxonomy method
    if reverse_ranks:
        # species, genus, family, order, class, phylum, superkingdom
        ranks = [rank for rank in NCBI.CANONICAL_RANKS]
    else:
        # superkingdom, phylum, class, order, family, genus, species
        ranks = [rank for rank in reversed(NCBI.CANONICAL_RANKS)]
    ranks.remove("root")
    starting_rank_index = ranks.index(starting_rank)
    ranks = ranks[starting_rank_index:]
    logger.debug(f"Using ranks: {', '.join(ranks)}")
    clustered_contigs = set()
    num_clusters = 0
    clusters = []
    for rank in ranks:
        # TODO: We should account for novel taxa here instead of removing 'unclassified'
        unclassified_filter = main[rank] != "unclassified"
        main_grouped_by_rank = main.loc[unclassified_filter].groupby(rank)
        taxa_counts = main_grouped_by_rank[rank].count()
        n_contigs_in_taxa = taxa_counts.sum()
        n_taxa = taxa_counts.index.nunique()
        logger.info(
            f"Examining {rank}: {n_taxa:,} unique taxa ({n_contigs_in_taxa:,} contigs)"
        )
        # Group contigs by rank and find best clusters within subset
        for rank_name_txt, dff in main_grouped_by_rank:
            if dff.empty:
                continue
            # Only cluster contigs that have not already been assigned a bin.
            # First filter with 'cluster' column
            rank_df = (
                dff.loc[dff["cluster"].isna()] if "cluster" in dff.columns else dff
            )
            # Second filter with previous clustering rounds' clustered contigs
            if clustered_contigs:
                rank_df = rank_df.loc[~rank_df.index.isin(clustered_contigs)]
            # After all of the filters, are there multiple contigs to cluster?
            if rank_df.empty:
                continue
            # Find best clusters
            logger.debug(
                f"Examining taxonomy: {rank} : {rank_name_txt} : {rank_df.shape}"
            )
            clusters_df = get_clusters(
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
            is_clustered = clusters_df["cluster"].notnull()
            clustered = clusters_df.loc[is_clustered]
            if clustered.empty:
                continue
            clustered_contigs.update({contig for contig in clustered.index})
            translation = {
                c: f"bin_{1+i+num_clusters:04d}"
                for i, c in enumerate(clustered.cluster.unique())
            }

            def rename_cluster(c):
                return translation[c]

            clustered.cluster = clustered.cluster.map(rename_cluster)
            num_clusters += clustered.cluster.nunique()
            clusters.append(clustered)

    if not clusters:
        raise BinningError("Failed to recover any clusters from dataset")
    # At this point we've went through all ranks and have clusters for each canonical-rank
    # We place these into one table and then...
    clustered_df = pd.concat(clusters, sort=True)
    # create a dataframe for any contigs *not* in the clustered dataframe
    unclustered_df = main.loc[~main.index.isin(clustered_df.index)]
    unclustered_df["cluster"] = pd.NA
    return pd.concat([clustered_df, unclustered_df], sort=True)


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
        help="Path to embedded k-mer frequencies table",
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
        "--domain",
        help="Kingdom to consider",
        choices=["bacteria", "archaea"],
        default="bacteria",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        default=False,
        help="log debug information",
    )
    args = parser.parse_args()
    # First merge all annotations
    kmers_df = pd.read_csv(args.kmers, sep="\t", index_col="contig")
    cov_df = pd.read_csv(args.coverages, sep="\t", index_col="contig")
    main_df = pd.merge(
        kmers_df,
        cov_df[["coverage"]],
        how="left",
        left_index=True,
        right_index=True,
    )
    gc_content_df = pd.read_csv(args.gc_content, sep="\t", index_col="contig")
    main_df = pd.merge(
        main_df,
        gc_content_df,
        how="left",
        left_index=True,
        right_index=True,
    )
    # Now read in markers for completeness and purity calculations
    markers_df = load_markers(args.markers)

    if args.taxonomy:
        taxa_df = pd.read_csv(args.taxonomy, sep="\t", index_col="contig")
        # Standardize naming (lower all rank names across taxonomy columns)
        for rank in NCBI.CANONICAL_RANKS:
            if rank not in taxa_df.columns:
                continue
            taxa_df[rank] = taxa_df[rank].map(lambda name: name.lower())
        # Now we are ready to filter by superkingdom
        taxa_df = taxa_df[taxa_df.superkingdom == args.domain]
        # merge taxa annotations with other annotations
        main_df = pd.merge(
            left=main_df,
            right=taxa_df,
            how="inner",
            left_index=True,
            right_index=True,
        )

    taxa_present = True if "taxid" in main_df else False
    main_df = main_df.convert_dtypes()
    # Sanity checks
    logger.debug(f"main_df shape: {main_df.shape}")

    main_out = binning(
        main=main_df,
        markers=markers_df,
        taxonomy=taxa_present,
        starting_rank=args.starting_rank,
        reverse_ranks=args.reverse_ranks,
        domain=args.domain,
        completeness=args.completeness,
        purity=args.purity,
        coverage_stddev=args.cov_stddev_limit,
        gc_content_stddev=args.gc_stddev_limit,
        method=args.clustering_method,
        verbose=args.verbose,
    )

    # Write out binning results with their respective binning metrics
    outcols = [
        "cluster",
        "completeness",
        "purity",
        "coverage_stddev",
        "gc_content_stddev",
    ]
    main_out[outcols].to_csv(
        args.output_binning, sep="\t", index=True, header=True, float_format="%.5f"
    )
    logger.info(f"Wrote binning results to {args.output_binning}")
    if args.output_main:
        # First after binning relevant assignments/metrics place contig physical annotations
        annotation_cols = ["coverage", "gc_content", "length"]
        outcols.extend(annotation_cols)
        if args.taxonomy:
            # Add in taxonomy columns if taxa are present
            # superkingdom, phylum, class, order, family, genus, species
            taxa_cols = [
                rank for rank in reversed(NCBI.CANONICAL_RANKS) if rank != "root"
            ]
            taxa_cols.append("taxid")
            # superkingdom, phylum, class, order, family, genus, species, taxid
            outcols.extend(taxa_cols)
        # Finally place kmer embeddings at end
        kmer_cols = [col for col in main_out.columns if "x_" in col]
        outcols.extend(kmer_cols)
        # Now write out table with columns sorted using outcols list
        main_out[outcols].to_csv(args.output_main, sep="\t", index=True, header=True)
        logger.info(f"Wrote main table to {args.output_main}")


if __name__ == "__main__":
    main()
