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

Cluster contigs recursively searching for bins with highest completeness and purity.
"""

import logging
import os
import shutil
import tempfile

import pandas as pd
import numpy as np

from Bio import SeqIO
from sklearn.cluster import DBSCAN
from hdbscan import HDBSCAN

from autometa.common.markers import Markers
from autometa.common import kmers

# TODO: This should be from autometa.common.kmers import Kmers
# So later we can simply/and more clearly do Kmers.load(kmers_fpath).embed(method)
from autometa.common.exceptions import TableFormatError
from autometa.taxonomy.ncbi import NCBI

pd.set_option("mode.chained_assignment", None)


logger = logging.getLogger(__name__)


def add_metrics(df, markers_df, domain="bacteria"):
    """Adds the completeness and purity metrics to each respective contig in df.

    Parameters
    ----------
    df : pd.DataFrame
        index='contig' cols=['x','y','coverage','cluster']
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
    clusters = dict(list(df.groupby("cluster")))
    metrics = {"purity": {}, "completeness": {}}
    for cluster, dff in clusters.items():
        contigs = dff.index.tolist()
        summed_markers = markers_df[markers_df.index.isin(contigs)].sum()
        is_present = summed_markers >= 1
        is_single_copy = summed_markers == 1
        nunique_markers = summed_markers[is_present].index.nunique()
        num_single_copy_markers = summed_markers[is_single_copy].index.nunique()
        completeness = nunique_markers / expected_number * 100
        # Protect from divide by zero
        if nunique_markers == 0:
            purity = pd.NA
        else:
            purity = num_single_copy_markers / nunique_markers * 100
        metrics["completeness"].update({cluster: completeness})
        metrics["purity"].update({cluster: purity})
    metrics_df = pd.DataFrame(metrics, index=clusters.keys())
    merged_df = pd.merge(df, metrics_df, left_on="cluster", right_index=True)
    return merged_df, metrics_df


def run_dbscan(df, eps, dropcols=["cluster", "purity", "completeness"]):
    """Run clustering on `df` at provided `eps`.

    Notes
    -----

        * documentation for sklearn `DBSCAN <https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html>`_
        * documentation for `HDBSCAN <https://hdbscan.readthedocs.io/en/latest/index.html>`_

    Parameters
    ----------
    df : pd.DataFrame
        Contigs with embedded k-mer frequencies as ['x','y'] columns and optionally 'coverage' column
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
    ValueError
        sets `usecols` and `dropcols` may not share elements

    """
    for col in dropcols:
        if col in df.columns:
            df.drop(columns=col, inplace=True)
    n_samples = df.shape[0]
    if n_samples == 1:
        clusters = pd.Series([pd.NA], index=df.index, name="cluster")
        return pd.merge(df, clusters, how="left", left_index=True, right_index=True)
    cols = ["x", "y"]
    if "coverage" in df.columns:
        cols.append("coverage")
    if np.any(df.isnull()):
        raise TableFormatError(
            f"df is missing {df.isnull().sum().sum()} kmer/coverage annotations"
        )
    X = df.loc[:, cols].to_numpy()
    clusterer = DBSCAN(eps=eps, min_samples=1, n_jobs=-1).fit(X)
    clusters = pd.Series(clusterer.labels_, index=df.index, name="cluster")
    return pd.merge(df, clusters, how="left", left_index=True, right_index=True)


def recursive_dbscan(
    table, markers_df, domain, completeness_cutoff, purity_cutoff, verbose=False,
):
    """Carry out DBSCAN, starting at eps=0.3 and continuing until there is just one
    group.

    Break conditions to speed up pipeline:
    Give up if we've got up to eps 1.3 and still no complete and pure clusters
    Often when you start at 0.3 there are zero complete and pure clusters, because
    the groups are too small. Later, some are found as the groups enlarge enough, but
    after it becomes zero again, it is a lost cause and we may as well stop. On the
    other hand, sometimes we never find any groups, so perhaps we should give up if
    by EPS 1.3 we never find any complete/pure groups.

    Parameters
    ----------
    table : pd.DataFrame
        Contigs with embedded k-mer frequencies as ['x','y','z'] columns and
        optionally 'coverage' column
    markers_df : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]
    domain : str
        Kingdom to determine metrics (the default is 'bacteria').
        choices=['bacteria','archaea']
    completeness_cutoff : float
        `completeness_cutoff` threshold to retain cluster (the default is 20.0).
    purity_cutoff : float
        `purity_cutoff` threshold to retain cluster (the default is 90.0).
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
            columns = ['x,'y','z','coverage','cluster','purity','completeness']
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
        median_completeness = metrics_df[completeness_filter & purity_filter][
            "completeness"
        ].median()
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
    complete_and_pure_df = best_df.loc[completeness_filter & purity_filter]
    unclustered_df = best_df.loc[~best_df.index.isin(complete_and_pure_df.index)]
    if verbose:
        logger.debug(f"Best completeness median: {best_median:4.2f}")
    logger.debug(
        f"clustered: {len(complete_and_pure_df)} unclustered: {len(unclustered_df)}"
    )
    return complete_and_pure_df, unclustered_df


def run_hdbscan(
    df,
    min_cluster_size,
    min_samples,
    cache_dir=None,
    dropcols=["cluster", "purity", "completeness"],
):
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
    for col in dropcols:
        if col in df.columns:
            df.drop(columns=col, inplace=True)
    n_samples = df.shape[0]
    if n_samples == 1:
        clusters = pd.Series([pd.NA], index=df.index, name="cluster")
        return pd.merge(df, clusters, how="left", left_index=True, right_index=True)
    cols = ["x", "y"]
    if "coverage" in df.columns:
        cols.append("coverage")
    if np.any(df.isnull()):
        raise TableFormatError(
            f"df is missing {df.isnull().sum().sum()} kmer/coverage annotations"
        )
    X = df.loc[:, cols].to_numpy()
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
    table, markers_df, domain, completeness_cutoff, purity_cutoff, verbose=False,
):
    """Recursively run HDBSCAN starting with defaults and iterating the min_samples
     and min_cluster_size until only 1 cluster is recovered.

    Parameters
    ----------
    table : pd.DataFrame
        Contigs with embedded k-mer frequencies as ['x','y','z'] columns and
        optionally 'coverage' column
    markers_df : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]
    domain : str
        Kingdom to determine metrics (the default is 'bacteria').
        choices=['bacteria','archaea']
    completeness_cutoff : float
        `completeness_cutoff` threshold to retain cluster (the default is 20.0).
    purity_cutoff : float
        `purity_cutoff` threshold to retain cluster (the default is 90.0).
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
            columns = ['x,'y','z','coverage','cluster','purity','completeness']
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
        median_completeness = metrics_df[completeness_filter & purity_filter][
            "completeness"
        ].median()
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
    complete_and_pure_df = best_df.loc[completeness_filter & purity_filter]
    unclustered_df = best_df.loc[~best_df.index.isin(complete_and_pure_df.index)]
    if verbose:
        logger.debug(f"Best completeness median: {best_median:4.2f}")
    logger.debug(
        f"clustered: {len(complete_and_pure_df)} unclustered: {len(unclustered_df)}"
    )
    return complete_and_pure_df, unclustered_df


def get_clusters(
    master_df,
    markers_df,
    domain="bacteria",
    completeness=20.0,
    purity=90.0,
    method="dbscan",
    verbose=False,
):
    """Find best clusters retained after applying `completeness` and `purity` filters.

    Parameters
    ----------
    master_df : pd.DataFrame
        index=contig,
        cols=['x','y','coverage']
    markers_df : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]
    domain : str
        Kingdom to determine metrics (the default is 'bacteria').
        choices=['bacteria','archaea'].
    completeness : float
        `completeness` threshold to retain cluster (the default is 20.).
    purity : float
        `purity` threshold to retain cluster (the default is 90.).
    method : str
        Description of parameter `method` (the default is 'dbscan').
        choices = ['dbscan','hdbscan']
    verbose : bool
        log stats for each recursive_dbscan clustering iteration

    Returns
    -------
    pd.DataFrame
        `master_df` with ['cluster','completeness','purity'] columns added
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
            master_df, markers_df, domain, completeness, purity, verbose=verbose,
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
        master_df = unclustered_df
    return pd.concat(clusters, sort=True)


def binning(
    master,
    markers,
    domain="bacteria",
    completeness=20.0,
    purity=90.0,
    taxonomy=True,
    starting_rank="superkingdom",
    method="dbscan",
    reverse_ranks=False,
    verbose=False,
):
    """Perform clustering of contigs by provided `method` and use metrics to
    filter clusters that should be retained via `completeness` and `purity`
    thresholds.

    Parameters
    ----------
    master : pd.DataFrame
        index=contig,
        cols=['x','y']
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
        Description of parameter `purity` (the default is 90.).
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
        master with ['cluster','completeness','purity'] columns added

    Raises
    -------
    TableFormatError
        No marker information is availble for contigs to be binned.
    """
    # First check needs to ensure we have markers available to check binning quality...
    if master.loc[master.index.isin(markers.index)].empty:
        raise TableFormatError(
            "No markers for contigs in table. Unable to assess binning quality"
        )

    logger.info(f"Using {method} clustering method")
    if not taxonomy:
        return get_clusters(
            master_df=master,
            markers_df=markers,
            domain=domain,
            completeness=completeness,
            purity=purity,
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
        unclassified_filter = master[rank] != "unclassified"
        master_grouped_by_rank = master.loc[unclassified_filter].groupby(rank)
        taxa_counts = master_grouped_by_rank[rank].count()
        n_contigs_in_taxa = taxa_counts.sum()
        n_taxa = taxa_counts.index.nunique()
        logger.info(
            f"Examining {rank}: {n_taxa:,} unique taxa ({n_contigs_in_taxa:,} contigs)"
        )
        # Group contigs by rank and find best clusters within subset
        for rank_name_txt, dff in master_grouped_by_rank:
            if dff.empty:
                continue
            # Only cluster contigs that have not already been assigned a bin.
            # First filter with 'cluster' column
            # Second filter with previous clustering rounds' clustered contigs
            rank_df = (
                dff.loc[dff["cluster"].isna()] if "cluster" in dff.columns else dff
            )
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
                master_df=rank_df,
                markers_df=markers,
                domain=domain,
                completeness=completeness,
                purity=purity,
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
    clustered_df = pd.concat(clusters, sort=True)
    unclustered_df = master.loc[~master.index.isin(clustered_df.index)]
    unclustered_df.loc[:, "cluster"] = pd.NA
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
    parser.add_argument("kmers", help="</path/to/kmers.tsv>")
    parser.add_argument("coverage", help="</path/to/coverages.tsv>")
    parser.add_argument("markers", help="</path/to/markers.tsv>")
    parser.add_argument("out", help="</path/to/output.tsv>")
    parser.add_argument("--embedded-kmers", help="</path/to/embedded_kmers.tsv>")
    parser.add_argument(
        "--embedding-method",
        help="Embedding method to use",
        choices=["bhsne", "sksne", "umap"],
        default="bhsne",
    )
    parser.add_argument(
        "--clustering-method",
        help="Clustering method to use",
        choices=["dbscan", "hdbscan"],
        default="dbscan",
    )
    parser.add_argument(
        "--completeness", help="<completeness cutoff>", default=20.0, type=float
    )
    parser.add_argument("--purity", help="<purity cutoff>", default=90.0, type=float)
    parser.add_argument("--taxonomy", help="</path/to/taxonomy.tsv>")
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
        " When --reverse-ranks given, contigs will be split in order of species, genus, family, order, class, phylum, superkingdom.",
    )
    parser.add_argument(
        "--domain",
        help="Kingdom to consider (archaea|bacteria)",
        choices=["bacteria", "archaea"],
        default="bacteria",
    )
    parser.add_argument(
        "--verbose", action="store_true", default=False, help="log debug information"
    )
    args = parser.parse_args()
    # Kmers.load().embed(method=args.embedded_kmers)
    kmers_df = kmers.embed(
        kmers=args.kmers, embedded=args.embedded_kmers, method=args.embedding_method
    )

    cov_df = pd.read_csv(args.coverage, sep="\t", index_col="contig")
    master_df = pd.merge(
        kmers_df, cov_df[["coverage"]], how="left", left_index=True, right_index=True,
    )

    markers_df = Markers.load(args.markers)
    markers_df = markers_df.convert_dtypes()
    # Taxonomy.load()
    if args.taxonomy:
        taxa_df = pd.read_csv(args.taxonomy, sep="\t", index_col="contig")
        taxa_df = taxa_df[taxa_df.superkingdom == args.domain]
        master_df = pd.merge(
            left=master_df,
            right=taxa_df,
            how="inner",
            left_index=True,
            right_index=True,
        )

    taxa_present = True if "taxid" in master_df else False
    master_df = master_df.convert_dtypes()
    logger.debug(f"master_df shape: {master_df.shape}")
    master_out = binning(
        master=master_df,
        markers=markers_df,
        taxonomy=taxa_present,
        starting_rank=args.starting_rank,
        reverse_ranks=args.reverse_ranks,
        domain=args.domain,
        completeness=args.completeness,
        purity=args.purity,
        method=args.clustering_method,
        verbose=args.verbose,
    )

    # Output table
    outcols = ["cluster", "completeness", "purity"]
    master_out[outcols].to_csv(args.out, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()
