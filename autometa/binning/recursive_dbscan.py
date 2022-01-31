#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

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
from numba import config


from autometa.common.markers import load as load_markers

from autometa.common.exceptions import TableFormatError, BinningError
from autometa.taxonomy.ncbi import NCBI
from autometa.binning.utilities import (
    write_results,
    read_annotations,
    add_metrics,
    filter_taxonomy,
    apply_binning_metrics_filter,
    reindex_bin_names,
    zero_pad_bin_names,
)

pd.set_option("mode.chained_assignment", None)

# See: https://numba.readthedocs.io/en/stable/user/threading-layer.html
# for more information on setting the threading layer. This is to prevent the warning
# Numba: Attempted to fork from a non-main thread, the TBB library may be in an invalid state in the child process
config.THREADING_LAYER = "safe"

logger = logging.getLogger(__name__)


def run_dbscan(
    df: pd.DataFrame,
    eps: float,
    n_jobs: int = -1,
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
    clusterer = DBSCAN(eps=eps, min_samples=1, n_jobs=n_jobs).fit(X)
    clusters = pd.Series(clusterer.labels_, index=df.index, name="cluster")
    return pd.merge(df, clusters, how="left", left_index=True, right_index=True)


def recursive_dbscan(
    table: pd.DataFrame,
    markers_df: pd.DataFrame,
    completeness_cutoff: float,
    purity_cutoff: float,
    coverage_stddev_cutoff: float,
    gc_content_stddev_cutoff: float,
    n_jobs: int = -1,
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
        binned_df = run_dbscan(df=table, eps=eps, n_jobs=n_jobs)
        df, metrics_df = add_metrics(df=binned_df, markers_df=markers_df)
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
    core_dist_n_jobs: int = -1,
) -> pd.DataFrame:
    """Run clustering on `df` at provided `min_cluster_size`.

    Notes
    -----

        * Reasoning for parameter: `cluster_selection_method <https://hdbscan.readthedocs.io/en/latest/parameter_selection.html#leaf-clustering>`_
        * Reasoning for parameters: `min_cluster_size and min_samples <https://hdbscan.readthedocs.io/en/latest/parameter_selection.html>`_
        * Documentation for `HDBSCAN <https://hdbscan.readthedocs.io/en/latest/index.html>`_

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

    core_dist_n_jobs: int
        Number of parallel jobs to run in core distance computations.
        For ``core_dist_n_jobs`` below -1, (n_cpus + 1 + core_dist_n_jobs) are used.

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
    # Ignore any errors raised from trying to drop previous 'cluster' in our df.
    df = df.drop(columns="cluster", errors="ignore")
    n_samples = df.shape[0]
    if n_samples == 1:
        clusters = pd.Series([pd.NA], index=df.index, name="cluster")
        return pd.merge(df, clusters, how="left", left_index=True, right_index=True)

    # NOTE: all of our kmer embedded columns should correspond from "x_1" to "x_{embedding_dimensions}"
    features_cols = [col for col in df.columns if "x_" in col or col == "coverage"]
    # Subset what will go into clusterer to only features (kmer and coverage information)
    features_df = df[features_cols]
    if np.any(features_df.isnull()):
        raise TableFormatError(
            f"df is missing {df.isnull().sum().sum()} kmer/coverage annotations"
        )
    # Fit and predict clusters
    clusters = HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method="leaf",
        allow_single_cluster=True,
        memory=cache_dir,
        core_dist_n_jobs=core_dist_n_jobs,
    ).fit_predict(features_df.to_numpy())
    clusters = pd.Series(clusters, index=df.index, name="cluster")
    # NOTE: HDBSCAN labels outliers with -1
    outlier_label = -1
    clusters = clusters.loc[clusters.ne(outlier_label)]
    return pd.merge(df, clusters, how="left", left_index=True, right_index=True)


def recursive_hdbscan(
    table: pd.DataFrame,
    markers_df: pd.DataFrame,
    completeness_cutoff: float,
    purity_cutoff: float,
    coverage_stddev_cutoff: float,
    gc_content_stddev_cutoff: float,
    n_jobs: int = -1,
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
            core_dist_n_jobs=n_jobs,
        )
        df, metrics_df = add_metrics(df=binned_df, markers_df=markers_df)
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
    completeness: float,
    purity: float,
    coverage_stddev: float,
    gc_content_stddev: float,
    method: str,
    n_jobs: int = -1,
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
    recursive_clusterers = {"dbscan": recursive_dbscan, "hdbscan": recursive_hdbscan}
    if method not in recursive_clusterers:
        raise ValueError(
            f"Method: {method} not a choice. choose between {recursive_clusterers.keys()}"
        )
    clusterer = recursive_clusterers[method]

    num_clusters = 0
    clusters = []
    # Continue until clusters are no longer being recovered
    # break when either clustered_df or unclustered_df is empty
    while True:
        clustered_df, unclustered_df = clusterer(
            table=main,
            markers_df=markers_df,
            completeness_cutoff=completeness,
            purity_cutoff=purity,
            coverage_stddev_cutoff=coverage_stddev,
            gc_content_stddev_cutoff=gc_content_stddev,
            verbose=verbose,
            n_jobs=n_jobs,
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
    df = pd.concat(clusters, axis="index", sort=True)
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


def taxon_guided_binning(
    main: pd.DataFrame,
    markers: pd.DataFrame,
    completeness: float = 20.0,
    purity: float = 95.0,
    coverage_stddev: float = 25.0,
    gc_content_stddev: float = 5.0,
    starting_rank: str = "superkingdom",
    method: str = "dbscan",
    reverse_ranks: bool = False,
    n_jobs: int = -1,
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

    completeness : float, optional
        Description of parameter `completeness` (the default is 20.).

    purity : float, optional
        purity threshold to retain cluster (the default is 95.0).
        e.g. cluster purity >= purity_cutoff

    coverage_stddev : float, optional
        cluster coverage threshold to retain cluster (the default is 25.0).

    gc_content_stddev : float, optional
        cluster GC content threshold to retain cluster (the default is 5.0).

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
    if reverse_ranks:
        # species, genus, family, order, class, phylum, superkingdom
        ranks = [rank for rank in NCBI.CANONICAL_RANKS if rank != "root"]
    else:
        # superkingdom, phylum, class, order, family, genus, species
        ranks = [rank for rank in reversed(NCBI.CANONICAL_RANKS) if rank != "root"]
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
                completeness=completeness,
                purity=purity,
                coverage_stddev=coverage_stddev,
                gc_content_stddev=gc_content_stddev,
                method=method,
                n_jobs=n_jobs,
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
        help="Path to embedded k-mers table",
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
        "--rank-filter",
        help="Taxonomy column canonical rank to subset by provided value of `--rank-name-filter`",
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
    parser.add_argument(
        "--cpus",
        default=-1,
        metavar="int",
        type=int,
        help="Number of cores to use by clustering method (default will try to use as many as are available)",
    )
    args = parser.parse_args()

    # First check if we are performing binning with taxonomic partitioning
    if args.taxonomy:
        main_df = read_annotations(
            [args.kmers, args.coverages, args.gc_content, args.taxonomy]
        )
        main_df = filter_taxonomy(
            df=main_df, rank=args.rank_filter, name=args.rank_name_filter
        )
    else:
        main_df = read_annotations([args.kmers, args.coverages, args.gc_content])

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

    if args.taxonomy:
        main_out = taxon_guided_binning(
            main=main_df,
            markers=markers_df,
            completeness=args.completeness,
            purity=args.purity,
            coverage_stddev=args.cov_stddev_limit,
            gc_content_stddev=args.gc_stddev_limit,
            method=args.clustering_method,
            starting_rank=args.starting_rank,
            reverse_ranks=args.reverse_ranks,
            n_jobs=args.cpus,
            verbose=args.verbose,
        )
    else:
        # Perform clustering w/o taxonomy
        main_out = get_clusters(
            main=main_df,
            markers_df=markers_df,
            completeness=args.completeness,
            purity=args.purity,
            coverage_stddev=args.cov_stddev_limit,
            gc_content_stddev=args.gc_stddev_limit,
            method=args.clustering_method,
            n_jobs=args.cpus,
            verbose=args.verbose,
        )

    write_results(
        results=main_out,
        binning_output=args.output_binning,
        full_output=args.output_main,
    )


if __name__ == "__main__":
    import sys

    # Using an http error status code...
    # From: https://kinsta.com/blog/http-status-codes/#200-status-codes
    # 204: “No Content.”
    # This code means that the server has successfully processed the request
    # but is not going to return any content.

    try:
        main()
    except (TableFormatError, BinningError) as err:
        logger.warn(err)
        sys.exit(204)
