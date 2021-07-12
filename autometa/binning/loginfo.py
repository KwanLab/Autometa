#!/usr/bin/env python

import argparse
import os
import re

import pandas as pd

from datetime import datetime
from typing import Dict, Union


def format_total_times(total_times: list, max_partition_size: str) -> pd.DataFrame:
    """Format total runtimes from timedelta objects to hours

    Parameters
    ----------
    total_times : list
        Runtime totals per algorithm during large-data-mode
    max_partition_size : str
        Partition size parameter retrieved from log file

    Returns
    -------
    pd.DataFrame
        Formatted dataframe of timedelta objects and hours per algorithm in large-data-mode run.
    """
    totals = pd.Series(total_times)
    totals["misc"] = totals.total_runtime - (
        totals.clustering + totals.kmer_count_normalization + totals.embedding
    )
    totals["total_hours"] = totals.total_runtime.total_seconds() / (60 * 60)
    totals["embedding_hours"] = totals.embedding.total_seconds() / (60 * 60)
    totals["clustering_hours"] = totals.clustering.total_seconds() / (60 * 60)
    totals[
        "kmer_count_normalization_hours"
    ] = totals.kmer_count_normalization.total_seconds() / (60 * 60)
    totals["misc_hours"] = totals.misc.total_seconds() / (60 * 60)
    totals.name = max_partition_size
    totals = totals.to_frame().transpose()
    totals.index.name = "max_partition_size"
    totals = totals.reset_index()
    return totals


def add_embedding_runtime_summary_info(
    embedding_df: pd.DataFrame, totals: pd.DataFrame
) -> pd.DataFrame:
    """Retrieve information about the embeddings that took the longest.

    Parameters
    ----------
    embedding_df : pd.DataFrame
        Embedding info retrieved from logfile
    totals : pd.DataFrame
        runtime totals summary table

    Returns
    -------
    pd.DataFrame
        runtime totals summary table updated with embedding info
    """

    ### embeddings for canonical ranks
    canonical_rank_embedding_df = embedding_df[
        embedding_df.taxon == embedding_df.canonical_rank
    ]
    totals["longest_canonical_rank_embedding_taxon"] = canonical_rank_embedding_df.loc[
        canonical_rank_embedding_df.duration.idxmax()
    ].taxon
    totals[
        "longest_canonical_rank_embedding_canonical_rank"
    ] = canonical_rank_embedding_df.loc[
        canonical_rank_embedding_df.duration.idxmax()
    ].canonical_rank
    totals[
        "longest_canonical_rank_embedding_num_contigs"
    ] = canonical_rank_embedding_df.loc[
        canonical_rank_embedding_df.duration.idxmax()
    ].num_data_points
    totals[
        "longest_canonical_rank_embedding_duration"
    ] = canonical_rank_embedding_df.loc[
        canonical_rank_embedding_df.duration.idxmax()
    ].duration
    ### embeddings for specific taxon
    taxon_embedding_df = embedding_df[embedding_df.taxon != embedding_df.canonical_rank]
    totals["longest_specific_taxon_embedding_taxon"] = taxon_embedding_df.loc[
        taxon_embedding_df.duration.idxmax()
    ].taxon
    totals["longest_specific_taxon_embedding_canonical_rank"] = taxon_embedding_df.loc[
        taxon_embedding_df.duration.idxmax()
    ].canonical_rank
    totals["longest_specific_taxon_embedding_num_contigs"] = taxon_embedding_df.loc[
        taxon_embedding_df.duration.idxmax()
    ].num_data_points
    totals["longest_specific_taxon_embedding_duration"] = taxon_embedding_df.loc[
        taxon_embedding_df.duration.idxmax()
    ].duration
    ### longest over all embeddings
    totals["longest_embedding_taxon"] = embedding_df.loc[
        embedding_df.duration.idxmax()
    ].taxon
    totals["longest_embedding_canonical_rank"] = embedding_df.loc[
        embedding_df.duration.idxmax()
    ].canonical_rank
    totals["longest_embedding_num_contigs"] = embedding_df.loc[
        embedding_df.duration.idxmax()
    ].num_data_points
    totals["longest_embedding_duration"] = embedding_df.loc[
        embedding_df.duration.idxmax()
    ].duration
    return totals


def add_clustering_runtime_summary_info(
    clustering_df: pd.DataFrame, totals: pd.DataFrame
) -> pd.DataFrame:
    """Retrieve information about the clustering that took the longest.

    Parameters
    ----------
    clustering_df : pd.DataFrame
        Clustering runtime info retrieved from logfile
    totals : pd.DataFrame
        runtime totals summary table

    Returns
    -------
    pd.DataFrame
        runtime totals summary table updated with clustering runtime info
    """
    ## Locate the taxon that took longest to cluster
    totals["longest_clustering_taxon"] = clustering_df.loc[
        clustering_df.duration.idxmax()
    ].taxon
    totals["longest_clustering_canonical_rank"] = clustering_df.loc[
        clustering_df.duration.idxmax()
    ].canonical_rank
    totals["longest_clustering_num_contigs"] = clustering_df.loc[
        clustering_df.duration.idxmax()
    ].num_data_points
    totals["longest_clustering_duration"] = clustering_df.loc[
        clustering_df.duration.idxmax()
    ].duration
    return totals


def get_runtimes(logfile: str) -> Dict[str, Union[pd.DataFrame, pd.Series, int, str]]:
    # TODO:
    # 1. Get total embedding time
    # 2. Num embeddings used that were already cached
    # 3. Get total k-mer count normalization time
    # 4. Get taxa skipped that had n.contigs above max_partition_size
    # 5. Get total binning time
    # Compile regular expressions
    max_partition_size_pattern = re.compile(r"Max\spartition\ssize\sset\sto:\s(\d+)")
    max_partition_size_write_pattern = re.compile(r"partitionSize\d+")
    begin_time_pattern = re.compile(
        r"\[(\d+/\d+/\d+\s\d+:\d+:\d+\s\S+)\sINFO\]\s\S+:\sSelected\sclustering\smethod:\s(\w+)"
    )
    skipped_taxa_pattern = re.compile(
        r"\[\d+/\d+/\d+\s\d+:\d+:\d+\s\S+\sDEBUG\]\s\S+:\s(\w+)\s>\smax_partition_size\s\((\d+,*\d*)>\w+"
    )
    cached_embedding_pattern = re.compile(r"Found\scached\sembeddings")
    norm_method_pattern = re.compile(
        r"\[(\d+/\d+/\d+\s\d+:\d+:\d+\s\S+)\sDEBUG\]\s\S+:\sTransforming\sk-mer\scounts\susing\s(\w+)"
    )
    pca_pattern = re.compile(
        r"\[(\d+/\d+/\d+\s\d+:\d+:\d+\s\S+)\sDEBUG\]\sautometa.common.kmers:\sPerforming\sdecomposition\swith\sPCA"
    )
    embed_begin_pattern = re.compile(
        r"\[(\d+/\d+/\d+\s\d+:\d+:\d+\s\S+)\sDEBUG\]\s\S+:\s(\w+):\s(\d+)\sdata\spoints\sand\s(\d+)\sdimensions"
    )
    embed_end_pattern = re.compile(
        r"\[(\d+/\d+/\d+\s\d+:\d+:\d+\s\S+)\sDEBUG\]\s\S+:\sCached\sembeddings\sto\s(\S+)"
    )
    cluster_begin_pattern = re.compile(
        r"\[(\d+/\d+/\d+\s\d+:\d+:\d+\s\S+)\sDEBUG\]\s\S+:\sExamining\staxonomy:\s(\w+)\s:\s(\w+)\s:\s\((\d+),\s\d+\)"
    )
    write_bins_pattern = re.compile(
        r"\[(\d+/\d+/\d+\s\d+:\d+:\d+\s\S+)\sINFO\]\s\S+:\sWrote\sbinning\sresults\sto"
    )
    # Now search log for patterns
    embeddings = []
    norm_transforms = []
    clusterings = []
    num_cached_embeddings = 0
    skipped_taxa = []
    total_times = {}
    cluster_begin_timestamp = None
    # These we only need to retrieve once from the log
    max_partition_size = None
    cluster_method = None
    embed_dims = None
    embed_method = None
    with open(logfile) as fh:
        for line in fh:

            # Get run begin timestamp and selected clustering method
            begin_time_match = begin_time_pattern.search(line)
            if begin_time_match:
                run_begin_timestamp = datetime.strptime(
                    begin_time_match.group(1), "%m/%d/%Y %I:%M:%S %p"
                )
                cluster_method = begin_time_match.group(2)

            # Get max partition size
            max_partition_size_match = max_partition_size_pattern.search(line)
            if max_partition_size_match and not max_partition_size:
                max_partition_size = max_partition_size_match.group(1)
                continue
            max_partition_size_write_match = max_partition_size_write_pattern.search(
                line
            )
            if max_partition_size_write_match and not max_partition_size:
                max_partition_size = max_partition_size_write_match.group().split(
                    "partitionSize"
                )[-1]

            # Find skipped embeddings
            cached_embedding_match = cached_embedding_pattern.search(line)
            if cached_embedding_match:
                num_cached_embeddings += 1
                continue

            # Find skipped taxa
            skipped_taxa_match = skipped_taxa_pattern.search(line)
            if skipped_taxa_match:
                skipped_taxon = skipped_taxa_match.group(1)
                n_contigs = int(skipped_taxa_match.group(2).replace(",", ""))
                skipped_taxa.append(
                    {"taxon": skipped_taxon, "num_data_points": n_contigs}
                )

            # Find beginning time of kmer normalization
            norm_method_match = norm_method_pattern.search(line)
            if norm_method_match:
                norm_method_begin_timestamp = datetime.strptime(
                    norm_method_match.group(1), "%m/%d/%Y %I:%M:%S %p"
                )
                norm_method = norm_method_match.group(2)
                norm_method_end_timestamp = None

            # Find end time of kmer normalization & begin time of PCA *OR* embedding
            pca_begin_match = pca_pattern.search(line)
            if pca_begin_match:
                norm_method_end_timestamp = datetime.strptime(
                    pca_begin_match.group(1), "%m/%d/%Y %I:%M:%S %p"
                )

            embed_begin_match = embed_begin_pattern.search(line)
            if embed_begin_match:
                embed_begin_timestamp = datetime.strptime(
                    embed_begin_match.group(1), "%m/%d/%Y %I:%M:%S %p"
                )
                # Determine kmer norm method duration
                ## perform this check to account for whether user performed PCA or not
                if not norm_method_end_timestamp:
                    norm_method_end_timestamp = embed_begin_timestamp
                norm_method_duration = (
                    norm_method_end_timestamp - norm_method_begin_timestamp
                )
                # Other metadata from embedding..
                if not embed_method:
                    embed_method = embed_begin_match.group(2)
                num_data_points = int(embed_begin_match.group(3).replace(",", ""))
                num_initial_dimensions = embed_begin_match.group(2)
                # Build norm_transforms table
                norm_transforms.append(
                    {
                        "begin_time": norm_method_begin_timestamp,
                        "end_time": norm_method_end_timestamp,
                        "duration": norm_method_duration,
                        "num_data_points": num_data_points,
                    }
                )
                if "kmer_count_normalization" not in total_times:
                    total_times["kmer_count_normalization"] = norm_method_duration
                else:
                    total_times["kmer_count_normalization"] += norm_method_duration

            # Find end time of kmer embedding
            embed_end_match = embed_end_pattern.search(line)
            if embed_end_match:
                embed_end_timestamp = datetime.strptime(
                    embed_end_match.group(1), "%m/%d/%Y %I:%M:%S %p"
                )
                embed_method_duration = embed_end_timestamp - embed_begin_timestamp
                embed_taxon = os.path.basename(embed_end_match.group(2)).split(
                    f".{norm_method}"
                )[0]
                embed_canonical_rank = os.path.basename(
                    os.path.dirname(embed_end_match.group(2))
                )
                if not embed_dims:
                    embed_dims = (
                        os.path.basename(embed_end_match.group(2))
                        .split(f"{embed_method}")[-1]
                        .split(".tsv")[0]
                    )
                embeddings.append(
                    {
                        "begin_time": embed_begin_timestamp,
                        "end_time": embed_end_timestamp,
                        "duration": embed_method_duration,
                        "taxon": embed_taxon,
                        "canonical_rank": embed_canonical_rank,
                        "num_data_points": num_data_points,  # this we use from embed_begin_match
                    }
                )
                if "embedding" not in total_times:
                    total_times["embedding"] = embed_method_duration
                else:
                    total_times["embedding"] += embed_method_duration

            # Find begin time of clustering
            cluster_begin_match = cluster_begin_pattern.search(line)
            if cluster_begin_match and not cluster_begin_timestamp:
                cluster_begin_timestamp = datetime.strptime(
                    cluster_begin_match.group(1), "%m/%d/%Y %I:%M:%S %p"
                )
                canonical_rank = cluster_begin_match.group(2)
                taxon = cluster_begin_match.group(3)
                num_taxa_contigs = int(cluster_begin_match.group(4).replace(",", ""))

            # end time of clustering is when we hit norm_method_begin_pattern ...
            if norm_method_match and cluster_begin_timestamp:
                cluster_end_timestamp = norm_method_begin_timestamp
                cluster_duration = cluster_end_timestamp - cluster_begin_timestamp
                clusterings.append(
                    {
                        "begin_time": cluster_begin_timestamp,
                        "end_time": cluster_end_timestamp,
                        "duration": cluster_duration,
                        "taxon": taxon,  # this we use from cluster_begin_match
                        "canonical_rank": canonical_rank,  # this we use from cluster_begin_match
                        "num_data_points": num_taxa_contigs,  # this we use from cluster_begin_match
                    }
                )
                cluster_begin_timestamp = None
                if "clustering" not in total_times:
                    total_times["clustering"] = cluster_duration
                else:
                    total_times["clustering"] += cluster_duration

            # or could reach end of taxa, so we write table signifying end of clustering
            write_bins_match = write_bins_pattern.search(line)
            if write_bins_match:
                run_end_timestamp = datetime.strptime(
                    write_bins_match.group(1), "%m/%d/%Y %I:%M:%S %p"
                )
                run_duration = run_end_timestamp - run_begin_timestamp
                total_times["total_runtime"] = run_duration
                # Add final clustering time in if the job was previously performing clustering...
                if cluster_begin_timestamp:
                    cluster_end_timestamp = run_end_timestamp
                    cluster_duration = cluster_end_timestamp - cluster_begin_timestamp
                    clusterings.append(
                        {
                            "begin_time": cluster_begin_timestamp,
                            "end_time": cluster_end_timestamp,
                            "duration": cluster_duration,
                            "taxon": taxon,  # this we use from cluster_begin_match
                            "canonical_rank": canonical_rank,  # this we use from cluster_begin_match
                            "num_data_points": num_taxa_contigs,  # this we use from cluster_begin_match
                        }
                    )
                    if "clustering" not in total_times:
                        total_times["clustering"] = cluster_duration
                    else:
                        total_times["clustering"] += cluster_duration
                break

    max_partition_size = "NA" if not max_partition_size else max_partition_size
    cluster_method = "NA" if not cluster_method else cluster_method

    norm_df = pd.DataFrame(norm_transforms)
    embedding_df = pd.DataFrame(embeddings)
    clustering_df = pd.DataFrame(clusterings)
    skipped_taxa_df = pd.DataFrame(skipped_taxa)

    num_total_embeddings = num_cached_embeddings + embedding_df.shape[0]

    totals = format_total_times(
        total_times=total_times, max_partition_size=max_partition_size
    )
    totals["embed_method"] = embed_method
    totals["embed_dims"] = embed_dims
    totals["num_total_embeddings"] = num_total_embeddings
    totals["num_cached_embeddings"] = num_cached_embeddings
    totals["cluster_method"] = cluster_method
    totals["num_taxa_skipped"] = skipped_taxa_df.shape[0]

    totals = add_embedding_runtime_summary_info(
        embedding_df=embedding_df, totals=totals
    )
    totals = add_clustering_runtime_summary_info(
        clustering_df=clustering_df, totals=totals
    )

    return {
        "embedding": embedding_df,
        "kmer_count_normalization": norm_df,
        "clustering": clustering_df,
        "skipped_taxa": skipped_taxa_df,
        "total": totals,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve clustering time stats from autometa.binning.recursive_dbscan err log"
    )
    parser.add_argument(
        "--log",
        help="Path to binning log file (If using slurm, this is typically stderr output path)",
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Directory to write runtime information tables",
        required=False,
        default=os.curdir,
    )
    parser.add_argument(
        "--prefix",
        help="Prefix to prepend to runtime information tables (Do not use a directory path as a prefix)",
        required=False,
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing log info table if it already exists",
        required=False,
        action="store_true",
        default=False,
    )
    args = parser.parse_args()
    # Retrieve runtimes from recursive_dbscan stderr logfile
    log_info = get_runtimes(args.log)
    # Make output directory if it does not exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # Write respecive runtime tables
    for info_type, info in log_info.items():
        if isinstance(info, pd.DataFrame):
            raw_filename = (
                f"{info_type}.tsv"
                if info_type == "skipped_taxa"
                else f"{info_type}_runtimes.tsv"
            )
            out_filename = (
                f"{args.prefix}{raw_filename}" if args.prefix else raw_filename
            )
            out_filepath = os.path.join(args.output, out_filename)
            if (
                os.path.exists(out_filepath)
                and os.path.getsize(out_filepath)
                and not args.overwrite
            ):
                print(
                    f"{out_filepath} already exists. User --overwrite to overwrite this table. Skipping..."
                )
            else:
                # In other scenarios (--overwrite provided or file does not exists) we write out table
                info.to_csv(out_filepath, sep="\t", index=False, header=True)
                print(f"Wrote {info_type} runtimes to {out_filepath}")
        else:
            print(f"{info_type}\t{info}")

    # Some misc summary methods:
    ## Retrieve total duration for individual algorithms.
    # embeddings.duration.sum()
    # norm_transforms.duration.sum()
    # clusterings.duration.sum()


if __name__ == "__main__":
    main()
