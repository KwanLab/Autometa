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

binning utilities script for autometa-binning

Script containing utility functions when performing autometa clustering/classification tasks.
"""


import logging

import pandas as pd

from typing import Iterable, Tuple

from autometa.taxonomy.ncbi import NCBI


logger = logging.getLogger(__name__)


def read_annotations(annotations: Iterable, how: str = "inner") -> pd.DataFrame:
    """Read in a list of contig annotations from filepaths and return all provided annotations in a single dataframe.

    Parameters
    ----------
    annotations : Iterable
        Filepaths of annotations. These should all contain a 'contig' column to be used as the index

    how : str, optional
        How to join the provided annotations. By default will take the 'inner' or intersection of all contigs from `annotations`.

    Returns
    -------
    pd.DataFrame
        index_col='contig', cols=[annotations, ...]
    """
    # Read in tables and concatenate along annotations axis
    df = pd.concat(
        [
            pd.read_csv(annotation, sep="\t", index_col="contig")
            for annotation in annotations
        ],
        axis="columns",
        join=how,
    ).convert_dtypes()
    logger.debug(f"merged annotations shape: {df.shape}")
    return df


def add_metrics(
    df: pd.DataFrame, markers_df: pd.DataFrame, domain: str = "bacteria"
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Adds cluster metrics to each respective contig in df.

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
        `df` with added cluster metrics columns=['completeness', 'purity', 'coverage_stddev', 'gc_content_stddev']
        pd.DataFrame(index=clusters,  cols=['completeness', 'purity', 'coverage_stddev', 'gc_content_stddev'])

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
    metrics = []
    if "cluster" in df.columns:
        for cluster, dff in df.groupby("cluster"):
            pfam_counts = markers_df[markers_df.index.isin(dff.index)].sum()
            is_present = pfam_counts >= 1
            is_single_copy = pfam_counts == 1
            nunique_markers = pfam_counts[is_present].count()
            num_single_copy_markers = pfam_counts[is_single_copy].count()
            completeness = nunique_markers / expected_number * 100
            # Protect from divide by zero
            if nunique_markers == 0:
                purity = pd.NA
            else:
                purity = num_single_copy_markers / nunique_markers * 100
            if dff.shape[0] <= 1:
                coverage_stddev = 0.0
                gc_content_stddev = 0.0
            else:
                coverage_stddev = dff.coverage.std()
                gc_content_stddev = dff.gc_content.std()
            metrics.append(
                {
                    "cluster": cluster,
                    "completeness": completeness,
                    "purity": purity,
                    "coverage_stddev": coverage_stddev,
                    "gc_content_stddev": gc_content_stddev,
                }
            )
    # Account for exceptions where clusters were not recovered
    if not metrics or "cluster" not in df.columns:
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

    else:
        metrics_df = pd.DataFrame(metrics).set_index("cluster")
        merged_df = pd.merge(df, metrics_df, left_on="cluster", right_index=True)
    return merged_df, metrics_df


def apply_binning_metrics_filter(
    df: pd.DataFrame,
    completeness_cutoff: float = 20.0,
    purity_cutoff: float = 95.00,
    coverage_stddev_cutoff: float = 25.0,
    gc_content_stddev_cutoff: float = 5.0,
) -> pd.DataFrame:
    """Filter `df` by provided cutoff values.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing binning metrics 'completeness', 'purity', 'coverage_stddev' and 'gc_content_stddev'

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

    Returns
    -------
    pd.DataFrame
        Cutoff filtered `df`

    Raises
    ------
    KeyError
        One of metrics to apply cutoff does not exist in the `df` columns

    """
    metrics = ["completeness", "purity", "coverage_stddev", "gc_content_stddev"]
    for metric in metrics:
        if metric not in df.columns:
            raise KeyError(f"{metric} not in `df` columns: {df.columns}")
    gt_completeness_cutoff = df["completeness"] >= completeness_cutoff
    gt_purity_cutoff = df["purity"] >= purity_cutoff
    lt_coverage_stddev_cutoff = df["coverage_stddev"] <= coverage_stddev_cutoff
    lt_gc_content_stddev_cutoff = df["gc_content_stddev"] <= gc_content_stddev_cutoff
    return df[
        gt_completeness_cutoff
        & gt_purity_cutoff
        & lt_coverage_stddev_cutoff
        & lt_gc_content_stddev_cutoff
    ]


def reindex_bin_names(
    df: pd.DataFrame, cluster_col: str = "cluster", initial_index: int = 0
) -> pd.DataFrame:
    """Re-index `cluster_col` using the provided `initial_index` as the initial index number
    then enumerating from this to the number of bins in `cluster_col` of `df`.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing `cluster_col`
    cluster_col : str, optional
        Cluster column to apply reindexing, by default "cluster"
    initial_index : int, optional
        Starting index number when reindexing, by default 0

    Note
    ----
    The bin names will start one number above the `initial_index` number provided. Therefore, the default
    behavior is to use 0 as the `initial_index` meaning the first bin name will be `bin_1`.

    Example
    -------
    .. code-block:: python

        >>>import pandas as pd
        >>>from autometa.binning.utilities import reindex_bin_names
        >>>df = pd.read_csv("binning.tsv", sep='\t', index_col='contig')
        >>>reindex_bin_names(df, cluster_col='cluster', initial_index=0)
                    cluster  completeness  purity  coverage_stddev  gc_content_stddev
        contig
        k141_1102126   bin_1     90.647482   100.0          1.20951           1.461658
        k141_110415    bin_1     90.647482   100.0          1.20951           1.461658
        k141_1210233   bin_1     90.647482   100.0          1.20951           1.461658
        k141_1227553   bin_1     90.647482   100.0          1.20951           1.461658
        k141_1227735   bin_1     90.647482   100.0          1.20951           1.461658
        ...              ...           ...     ...              ...                ...
        k141_999969      NaN           NaN     NaN              NaN                NaN
        k141_99997       NaN           NaN     NaN              NaN                NaN
        k141_999982      NaN           NaN     NaN              NaN                NaN
        k141_999984      NaN           NaN     NaN              NaN                NaN
        k141_999987      NaN           NaN     NaN              NaN                NaN

    Returns
    -------
    pd.DataFrame
        [description]
    """
    # Create mapping from current cluster_col round name to padded name
    reindexed_bin_names = {
        bin_name: f"bin_{1+bin_num+initial_index}"
        for bin_num, bin_name in enumerate(df[cluster_col].dropna().unique())
    }
    # Define function to apply mapping
    def reindex_cluster(bin_name):
        return reindexed_bin_names[bin_name]

    # Apply mapping
    df[cluster_col] = df[cluster_col].dropna().map(reindex_cluster)
    return df


def zero_pad_bin_names(df: pd.DataFrame, cluster_col: str = "cluster") -> pd.DataFrame:
    """Apply zero padding to `cluster_col` using the length of digit corresponding to
    the number of unique clusters in `cluster_col` in the `df`.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe containing `cluster_col`
    cluster_col : str, optional
        Cluster column to apply zero padding, by default "cluster"

    Returns
    -------
    pd.DataFrame
        Dataframe with `cluster_col` zero padded to the length of the number of clusters
    """
    left_pad_zeros = len(str(df[cluster_col].dropna().nunique()))
    # Create mapping from current cluster_col round name to padded name
    zero_padded_bin_names = {
        # NOTE: We add 1 to the bin_num because python is 0-indexed and we want to start at bin num of 1
        bin_name: (
            f"bin_" + str(bin_num + 1).zfill(left_pad_zeros)
            # We need to account for the 'unclustered' contigs
            if not isinstance(bin_name, type(pd.NA)) or "bin_" in bin_name
            else bin_name
        )
        for bin_num, bin_name in enumerate(df[cluster_col].dropna().unique())
    }
    # Define function to apply mapping
    def zero_pad_cluster(cluster_name):
        return zero_padded_bin_names[cluster_name]

    # Apply bin name mapping (only to contigs assigned a cluster)
    df[cluster_col] = df[cluster_col].dropna().map(zero_pad_cluster)
    return df


def write_results(
    results: pd.DataFrame, binning_output: str, full_output: str = None
) -> None:
    """Write out binning results with their respective binning metrics

    Parameters
    ----------
    results : pd.DataFrame
        Binning results contigs dataframe consisting of "cluster" assignments with their respective metrics and annotations

    binning_output : str
        Filepath to write binning results

    full_output : str, optional
        If provided, will write assignments, metrics and annotations together into `full_output` (filepath)

    Returns
    -------
    NoneType

    """
    outcols = [
        "cluster",
        "completeness",
        "purity",
        "coverage_stddev",
        "gc_content_stddev",
    ]
    results[outcols].to_csv(binning_output, sep="\t", index=True, header=True)
    logger.info(f"Wrote binning results to {binning_output}")
    if full_output:
        # First after binning relevant assignments/metrics place contig physical annotations
        annotation_cols = ["coverage", "gc_content", "length"]
        outcols.extend(annotation_cols)
        # Add in taxonomy columns if taxa are present
        # superkingdom, phylum, class, order, family, genus, species
        taxa_cols = [rank for rank in reversed(NCBI.CANONICAL_RANKS) if rank != "root"]
        taxa_cols.append("taxid")
        # superkingdom, phylum, class, order, family, genus, species, taxid
        for taxa_col in taxa_cols:
            if taxa_col in results.columns:
                outcols.append(taxa_col)
        # Finally place kmer embeddings at end
        kmer_cols = [col for col in results.columns if "x_" in col]
        outcols.extend(kmer_cols)
        # Now write out table with columns sorted using outcols list
        results[outcols].to_csv(full_output, sep="\t", index=True, header=True)
        logger.info(f"Wrote main table to {full_output}")


if __name__ == "__main__":
    import sys

    print("Utilities script for autometa binning. This should not be run directly!")
    sys.exit(0)
