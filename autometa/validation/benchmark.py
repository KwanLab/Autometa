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

Autometa binning benchmarking.

Script to benchmark Autometa clustering results using clustering evaluation metrics.
# Setting seed? See: https://stackoverflow.com/a/5837352/13118765
"""


import logging

import os
from typing import Dict
import pandas as pd

from sklearn import metrics

logger = logging.getLogger(__name__)


def evaluate_clustering_performance(
    predictions: str, reference: str
) -> Dict[str, float]:
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

    reference : str
        Path to known reference community assignments to compare against Autometa binning results.

    Returns
    -------
    Dict[str, float]
        computed clustering evaluation metrics keyed by their metric

    Raises
    -------
    TableFormatError
        The provided reference community and predictions do not match!

    """
    pred_df = pd.read_csv(
        predictions, sep="\t", index_col="contig", usecols=["contig", "cluster"]
    )
    ref_df = pd.read_csv(reference, sep="\t", index_col="contig")
    # Modification here retrieves only contigs in binning dataset. This is not giving us the "complete" score.
    ref_df = ref_df[ref_df.index.isin(pred_df.index)]
    unclustered_idx = pred_df[pred_df.cluster == "unclustered"].index
    pred_df.loc[unclustered_idx, "cluster"] = pd.NA
    pred_df.cluster = pred_df.cluster.astype("category")
    # Remove any contigs that were assigned as misassembled
    ref_df = ref_df[ref_df.reference_genome != "misassembled"]
    ref_df.reference_genome = ref_df.reference_genome.astype("category")
    mdf = pd.merge(pred_df, ref_df, how="left", left_index=True, right_index=True)
    labels_true = mdf.reference_genome.cat.codes.values
    labels_pred = mdf.cluster.cat.codes.values
    ari = metrics.adjusted_rand_score(labels_true=labels_true, labels_pred=labels_pred)
    homogeneity_score = metrics.homogeneity_score(labels_true, labels_pred)
    completeness_score = metrics.completeness_score(labels_true, labels_pred)
    v_measure = metrics.v_measure_score(labels_true, labels_pred)
    ami = metrics.adjusted_mutual_info_score(
        labels_true=labels_true, labels_pred=labels_pred
    )
    gnmi = metrics.normalized_mutual_info_score(
        labels_true=labels_true,
        labels_pred=labels_pred,
        average_method="geometric",
    )
    fowlkes_mallows_score = metrics.fowlkes_mallows_score(
        labels_true=labels_true, labels_pred=labels_pred
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
    # Note: If you do not have any defaults corresponding to your parameters,
    # you may remove the formatter class: ArgumentDefaultsHelpFormatter
    # to reduce help text verbosity.
    parser = argparse.ArgumentParser(
        description="Concise description of script.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--binning",
        help="Path to Autometa binning output",
        required=True,
        nargs="*",
    )
    parser.add_argument(
        "--assignments",
        help="Path to community reference assignments",
        required=True,
        nargs="*",
    )
    parser.add_argument(
        "--output",
        help="Path to write clustering evaluation metrics",
        required=False,
        default="clustering_evaluation_metrics.tsv",
    )
    parser.add_argument(
        "--stacked",
        help="Path to write stacked clustering evaluation metrics",
        required=False,
        default="clustering_evaluation_metrics.stack.tsv",
    )
    args = parser.parse_args()
    metrics = []
    for binning, assignments in zip(args.binning, args.assignments):
        # Naming of dataset currently follows simulated_communities test data directory structure
        # Located here: https://drive.google.com/open?id=1JFjVb-pfQTv4GXqvqRuTOZTfKdT0MwhN
        dataset = os.path.basename(os.path.realpath(os.path.dirname(binning)))
        assert dataset in assignments, "Assignments do not match predictions!"
        metric = evaluate_clustering_performance(
            predictions=binning, reference=assignments
        )
        # dataset = "/path/to/dataset/binning.tsv.gz"
        metric.update({"dataset": dataset})
        metrics.append(metric)

    df = pd.DataFrame(metrics).set_index("dataset")
    # e.g. index = 78.125Mbp
    # TODO: Need to account for synthetic communities... i.e. MIX-51, MIX-51-equal, etc.
    df.sort_index(
        key=lambda x: x.str.replace("Mbp", "").astype(float),
        ascending=True,
        inplace=True,
    )
    df.to_csv(args.output, sep="\t", index=True, header=True)
    logger.info(f"Wrote {df.shape[0]} datasets metrics to {args.output}")
    # Write out stacked dataframe for visualization with `plot-cluster-evaluation-metrics.R`
    dff = df.stack()
    dff.index.name = ("dataset", "metric")
    dff.name = "score"
    dff = dff.to_frame().reset_index(level=1).rename(columns={"level_1": "metric"})
    dff.to_csv(args.stacked, sep="\t", index=True, header=True)
    logger.info(
        f"Wrote {dff.index.nunique()} datasets (stacked) metrics to {args.stacked}"
    )


if __name__ == "__main__":
    main()