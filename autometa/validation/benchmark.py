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

Autometa clustering evaluation benchmarking.

Script to benchmark Autometa clustering results using clustering evaluation metrics.
# Setting seed? See: https://stackoverflow.com/a/5837352/13118765
"""


import logging

import os
from typing import Dict, NamedTuple, Union
import pandas as pd
from collections import namedtuple

from sklearn import metrics

logger = logging.getLogger(__name__)

Labels = namedtuple("Labels", ["true", "pred"])


def get_categorical_labels(
    predictions: str, reference: Union[str, pd.DataFrame]
) -> NamedTuple:
    pred_df = pd.read_csv(
        predictions, sep="\t", index_col="contig", usecols=["contig", "cluster"]
    )
    if not isinstance(reference, pd.DataFrame) and isinstance(reference, str):
        ref_df = pd.read_csv(reference, sep="\t", index_col="contig")
    elif not isinstance(reference, pd.DataFrame) and not isinstance(reference, str):
        raise ValueError(f"reference is an invalid argument type: {type(reference)}")
    else:
        ref_df = reference
    # Modification here retrieves only contigs in binning dataset. This is not giving us the "complete" score.
    ref_df = ref_df[ref_df.index.isin(pred_df.index)]
    # Remove any contigs that were assigned as misassembled
    ref_df = ref_df[ref_df.reference_genome != "misassembled"]
    # Set reference_genome as category type so we can easily retrieve categorical labels
    # NOTE: These do not need to be a one-to-one correspondence to the prediction dataframe's categorical labels
    # (The exact label value does not matter for these metrics as we are only looking at the groupings [not classification])
    ref_df.reference_genome = ref_df.reference_genome.astype("category")
    # Assign "unclustered" to NA and convert 'cluster' column to categorical type
    unclustered_idx = pred_df[pred_df.cluster == "unclustered"].index
    pred_df.loc[unclustered_idx, "cluster"] = pd.NA
    pred_df.cluster = pred_df.cluster.astype("category")
    # Merge reference_assignments and predictions
    master_df = pd.merge(pred_df, ref_df, how="left", left_index=True, right_index=True)
    if master_df.empty:
        raise ValueError(
            "The provided reference community and predictions do not match!"
        )
    # Retrieve categorical values for each set of labels (truth and predicted)
    labels_true = master_df.reference_genome.cat.codes.values
    labels_pred = master_df.cluster.cat.codes.values
    return Labels(true=labels_true, pred=labels_pred)


def evaluate_clustering_performance(labels: NamedTuple) -> Dict[str, float]:
    """Calculate various clustering performance metrics listed below.

    Note
    ----

    Some of these clustering performance evaluation metrics adjust for chance. This is dicussed in more detail in
    the sklearn user guide. - `Adjustment for chance in clustering performance evaluation <https://scikit-learn.org/stable/auto_examples/cluster/plot_adjusted_for_chance_measures.html#sphx-glr-auto-examples-cluster-plot-adjusted-for-chance-measures-py>`_

    This analysis suggests the most robust and reliable metrics to use as the number of clusters
    increases are adjusted rand index and adjusted mutual info score.

    Metrics
    -------

    * `adjusted mutual info score <https://scikit-learn.org/stable/modules/clustering.html#mutual-information-based-scores>`_
    * `geometric normalized mutual info score <https://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html#sklearn-metrics-normalized-mutual-info-score>`_
    * `adjusted rand index <https://scikit-learn.org/stable/modules/clustering.html#rand-index>`_
    * `homogeneity <https://scikit-learn.org/stable/modules/generated/sklearn.metrics.homogeneity_score.html#sklearn.metrics.homogeneity_score>`_
    * `completeness <https://scikit-learn.org/stable/modules/generated/sklearn.metrics.completeness_score.html#sklearn.metrics.completeness_score>`_
    * `V-measure <https://scikit-learn.org/stable/modules/clustering.html#homogeneity-completeness-and-v-measure>`_
    * `fowlkes-mallows score <https://scikit-learn.org/stable/modules/clustering.html#fowlkes-mallows-scores>`_


    Parameters
    ----------
    predictions : str
        Path to Autometa binning results.

    reference : str|pd.DataFrame
        Path to or dataframe of known reference community assignments to compare against Autometa binning results.

    Returns
    -------
    Dict[str, float]
        computed clustering evaluation metrics keyed by their metric

    Raises
    -------
    ValueError
        The input arguments are not the correct type (pd.DataFrame or str)
    ValueError
        The provided reference community and predictions do not match!

    """
    # Calculate cluster variety of cluster evaluation metrics with labels
    ari = metrics.adjusted_rand_score(labels_true=labels.true, labels_pred=labels.pred)
    homogeneity_score = metrics.homogeneity_score(labels.true, labels.pred)
    completeness_score = metrics.completeness_score(labels.true, labels.pred)
    v_measure = metrics.v_measure_score(labels.true, labels.pred)
    ami = metrics.adjusted_mutual_info_score(
        labels_true=labels.true, labels_pred=labels.pred
    )
    gnmi = metrics.normalized_mutual_info_score(
        labels_true=labels.true,
        labels_pred=labels.pred,
        average_method="geometric",
    )
    fowlkes_mallows_score = metrics.fowlkes_mallows_score(
        labels_true=labels.true, labels_pred=labels.pred
    )
    return {
        "adjusted mutual info score": ami,
        "geometric normalized mutual info score": gnmi,
        "adjusted rand score": ari,
        "homogeneity score": homogeneity_score,
        "completeness score": completeness_score,
        "V-measure": v_measure,
        "fowlkes-mallows score": fowlkes_mallows_score,
    }


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Benchmark clustering against reference assignments for the provided simulated/synthetic community.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--binning",
        help="Path to Autometa binning output (May specify multiple if they all correspond to the same `--assignments` community ",
        metavar="filepath",
        nargs="*",
    )
    parser.add_argument(
        "--assignments",
        help="Path to community reference assignments",
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--output-wide",
        help="Path to write clustering evaluation metrics (each metric receives its own column)",
        metavar="filepath",
        required=False,
        default="clustering_benchmarks.tsv.gz",
    )
    parser.add_argument(
        "--output-long",
        help="Path to write clustering evaluation metrics (metrics are stacked into one 'metric' column)",
        metavar="filepath",
        required=False,
        default=None,
    )
    args = parser.parse_args()
    all_metrics = []
    for binning in args.binning:
        labels = get_categorical_labels(predictions=binning, reference=args.assignments)
        metrics = evaluate_clustering_performance(labels)
        metrics.update({"dataset": os.path.basename(binning)})
        all_metrics.append(metrics)

    df = pd.DataFrame(all_metrics).set_index("dataset")

    df.to_csv(args.output_wide, sep="\t", index=True, header=True)
    logger.info(f"Wrote {df.shape[0]} datasets metrics to {args.output_wide}")
    if args.output_long:
        # Write out stacked dataframe for visualization with `plot-cluster-evaluation-metrics.R`
        dff = df.stack()
        dff.index.name = ("dataset", "metric")
        dff.name = "score"
        dff = dff.to_frame().reset_index(level=1).rename(columns={"level_1": "metric"})
        dff.to_csv(args.output_long, sep="\t", index=True, header=True)
        logger.info(
            f"Wrote {dff.index.nunique()} datasets (stacked) metrics to {args.output_long}"
        )


if __name__ == "__main__":
    main()