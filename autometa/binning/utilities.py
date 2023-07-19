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

from autometa.taxonomy.database import TaxonomyDatabase


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
    logger.debug(f"Reading/merging {len(annotations):,} contig annotation files")
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
    for canonical_rank in TaxonomyDatabase.CANONICAL_RANKS:
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


def add_metrics_dynamic(
    df: pd.DataFrame,
    markers_df: pd.DataFrame,
    tax_marker_sets_type: str = 'single',
    contig_at_tax_df: pd.DataFrame,
    tax_marker_set_dict=dict,
    tax_assignment_dict=dict,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Adds cluster metrics to each respective contig in df using metrics calculations to accommodate for dynamic gene marker sets

     * runs LCA to each cluster to determine which marker genes to use *
     * note to shane - make sure this is taking in the tax_id column from an earlier step *




    Parameters
     ----------
     df : pd.DataFrame
         index='contig' cols=['coverage','gc_content','cluster','x_1','x_2',...,'x_n']

     markers_df : pd.DataFrame
         wide format, i.e. index=contig cols=[marker,marker,...]

     tax_marker_sets_type: str
         will calculate metrics using either single or multi count marker genes.  ["single", "multi"]

     contig_at_tax_df : pd.DataFrame
         contig_at_rank_df, i.e. index=contig cols=ranks

     tax_marker_set_dict: dict
         {tax_id: [gene_set_to_use]}

    tax_assignment_dict: dict
        {bin: tax_id}
         
     Returns
     -------
     2-tuple
         `df` with added cluster metrics columns=['completeness', 'purity', 'coverage_stddev', 'gc_content_stddev']
         pd.DataFrame(index=clusters,  cols=['completeness', 'purity', 'coverage_stddev', 'gc_content_stddev'])

    """

    metrics = [
        "completeness",
        "purity",
        "coverage_stddev",
        "gc_content_stddev",
        "present_marker_count",
        "single_copy_marker_count",
    ]

    if "cluster" not in df.columns:
        cluster_metrics_df = pd.DataFrame(
            data={metric: pd.NA for metric in metrics}, index=df.index
        )
        # Remove previous metrics to avoid creating metrics with suffixes
        contig_metrics_df = (
            df.copy().convert_dtypes().drop(columns=metrics, errors="ignore")
        )
        contig_metrics_df[metrics] = pd.NA

    else:

        # Remove previous metrics to avoid creating metrics with suffixes
        df = df.drop(columns=metrics, errors="ignore")
        # join cluster and marker data -> group contigs by cluster
        main_grouped_by_cluster = df.join(markers_df, how="outer").groupby("cluster")


        # shane make sure all markers are being used.
        # double check tax_bin assignment


        if tax_marker_sets_type == "single":
            # marker counts per cluster

            cluster_marker_counts = main_grouped_by_cluster[markers_df.columns].sum()
            # ie index= cluster 1, col= markers


            # number of single count markers per cluster
            single_copy_marker_count = cluster_marker_counts.eq(1).sum(axis=1)

            # this must be for only those in gene set
            # present_marker_count = cluster_marker_counts.ge(1).sum(axis=1)
            # [df[test_dict[key]].sum(axis=1).loc[key] for key in test_dict.keys()]
           
            present_marker_count = pd.Series()
            reference_marker_counts= pd.Series()
            single_copy_marker_count = pd.Series()
            
            # assumes df_tax_assignment = dict
            # gene_set_dict = dict

            for cluster in cluster_marker_counts.index:
                # get tax for this contig
                tax = tax_assignment_dict[cluster]
                # find which marker genes to look for
                marker_gene_set = tax_marker_set_dict[tax]
                # get counts for this contig for this gene set
                counts_in_this_cluster = cluster_marker_counts[marker_gene_set].loc[cluster]
                # number of marker genes present
                present_marker_count[cluster] = counts_in_this_cluster.ge(1)
                # number of single copy marker genes present
                single_copy_marker_count[cluster] = counts_in_this_cluster.eq(1)
                # number of marker genes in reference set
                reference_marker_counts[cluster] = len(marker_gene_set)


                

                ## cluster_marker_counts[marker_gene_set].loc[contig].sum(axis=1)

            # need to get lookup for this
            # pd.Series(test_list).map(test_dict).to_list()
            # [len(x) for x in pd.Series(test_list).map(test_dict).values]

            reference_marker_counts = [
                len(marker_list)
                for marker_list in cluster_marker_counts.index.map(
                    tax_gene_set_dict
                ).to_list()
            ]

            completeness = present_marker_count / reference_marker_counts * 100
            purity = single_copy_marker_count / present_marker_count * 100

            coverage_stddev = main_grouped_by_cluster.coverage.std()
            gc_content_stddev = main_grouped_by_cluster.gc_content.std()

        else:
            cluster_marker_counts = main_grouped_by_cluster[
                main_grouped_by_cluster[markers_df.columns] >= 1
            ].counts()

        # NOTE: df.ge(...) and df.eq(...) operators return boolean pd.DataFrame
        # count single-copy
        

        cluster_metrics_df = pd.DataFrame(
            {
                "completeness": completeness,
                "purity": purity,
                "coverage_stddev": coverage_stddev,
                "gc_content_stddev": gc_content_stddev,
                "present_marker_count": present_marker_count,
                "single_copy_marker_count": single_copy_marker_count,
            }
        )


        # TODO: add redundant marker count to mag summary
        # ... or 1, 2, 3, 4, 5+ count cols similar to CheckM output table
        # redundant_marker_count = cluster_marker_counts.gt(1).sum(axis=1)
        # calculate completeness and purity and std. dev. metrics

        contig_metrics_df = pd.merge(
            df, cluster_metrics_df, how="left", left_on="cluster", right_index=True
        )
    return contig_metrics_df, cluster_metrics_df


def add_metrics(
    df: pd.DataFrame, markers_df: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Adds cluster metrics to each respective contig in df.

    Metrics
    -------

        * :math:`completeness = \frac{markers_{cluster}}{markers_{ref}} * 100`
        * :math:`purity \% = \frac{markers_{single-copy}}{markers_{cluster}} * 100`
        * :math:`\mu_{coverage}  = \frac{1}{N}\sum_{i=1}^{N}\left(x_{i}-\mu\right)^{2}`
        * :math:`\mu_{GC\%}  = \frac{1}{N}\sum_{i=1}^{N}\left(x_{i}-\mu\right)^{2}`

    Parameters
    ----------
    df : pd.DataFrame
        index='contig' cols=['coverage','gc_content','cluster','x_1','x_2',...,'x_n']

    markers_df : pd.DataFrame
        wide format, i.e. index=contig cols=[marker,marker,...]

    tax_aware : Bool
        will use different metrics calculations to accommodate for tax aware marker sets (default: False)
    Returns
    -------
    2-tuple
        `df` with added cluster metrics columns=['completeness', 'purity', 'coverage_stddev', 'gc_content_stddev']
        pd.DataFrame(index=clusters,  cols=['completeness', 'purity', 'coverage_stddev', 'gc_content_stddev'])

    """
    metrics = [
        "completeness",
        "purity",
        "coverage_stddev",
        "gc_content_stddev",
        "present_marker_count",
        "single_copy_marker_count",
    ]
    reference_markers_count = markers_df.shape[1]
    # Account for exceptions where clusters were not recovered
    if "cluster" not in df.columns:
        cluster_metrics_df = pd.DataFrame(
            data={metric: pd.NA for metric in metrics}, index=df.index
        )
        # Remove previous metrics to avoid creating metrics with suffixes
        contig_metrics_df = (
            df.copy().convert_dtypes().drop(columns=metrics, errors="ignore")
        )
        contig_metrics_df[metrics] = pd.NA
    else:
        # Remove previous metrics to avoid creating metrics with suffixes
        df = df.drop(columns=metrics, errors="ignore")
        # join cluster and marker data -> group contigs by cluster
        main_grouped_by_cluster = df.join(markers_df, how="outer").groupby("cluster")
        # sum cluster markers
        cluster_marker_counts = main_grouped_by_cluster[markers_df.columns].sum()
        # NOTE: df.ge(...) and df.eq(...) operators return boolean pd.DataFrame
        # count single-copy
        single_copy_marker_count = cluster_marker_counts.eq(1).sum(axis=1)
        present_marker_count = cluster_marker_counts.ge(1).sum(axis=1)
        # TODO: add redundant marker count to mag summary
        # ... or 1, 2, 3, 4, 5+ count cols similar to CheckM output table
        # redundant_marker_count = cluster_marker_counts.gt(1).sum(axis=1)
        # calculate completeness and purity and std. dev. metrics
        completeness = present_marker_count / reference_markers_count * 100
        purity = single_copy_marker_count / present_marker_count * 100
        coverage_stddev = main_grouped_by_cluster.coverage.std()
        gc_content_stddev = main_grouped_by_cluster.gc_content.std()
        # merge metrics with given dataframe
        cluster_metrics_df = pd.DataFrame(
            {
                "completeness": completeness,
                "purity": purity,
                "coverage_stddev": coverage_stddev,
                "gc_content_stddev": gc_content_stddev,
                "present_marker_count": present_marker_count,
                "single_copy_marker_count": single_copy_marker_count,
            }
        )
        contig_metrics_df = pd.merge(
            df, cluster_metrics_df, how="left", left_on="cluster", right_index=True
        )
    return contig_metrics_df, cluster_metrics_df


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
    # df[metric].ge(cutoff) equivalent operator: metric >= cutoff
    # df[metric].le(cutoff) equivalent operator: metric <= cutoff
    return df[
        df.completeness.ge(completeness_cutoff)
        & df.purity.ge(purity_cutoff)
        & df.coverage_stddev.le(coverage_stddev_cutoff)
        & df.gc_content_stddev.le(gc_content_stddev_cutoff)
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
        DataFrame of re-indexed bins in `cluster_col` starting at `initial_index` + 1
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
        taxa_cols = [
            rank
            for rank in reversed(TaxonomyDatabase.CANONICAL_RANKS)
            if rank != "root"
        ]
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
