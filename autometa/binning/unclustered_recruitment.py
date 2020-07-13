#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth
Uppal, Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify it under the
terms of the GNU Affero General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version.

Autometa is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License along
with Autometa. If not, see <http://www.gnu.org/licenses/>. COPYRIGHT


Recruit unclustered contigs using metagenome annotations and binning results.

## Pseudocode for unclustered recruitment algorithm:

1. Retrieve all features from all contigs and load markers
2. Split contigs by clustered/unclustered
3. Subset clustered contigs by contigs containing markers
4. subset features by clustered (subset) and unclustered contigs
5. split clustered features into training and test set for estimator performance measure
    - split clustered contig features into a training and test set (50% random subset via `train_test_split`)
6. Get trained classifier's predictions on unclustered contigs' features
7. Jump back to 5 for provided number of classifications
8. Filter predictions:
    - Keep predictions if prediction count >= confidence threshold
        confidence threshold = num. consistent classifications / num. classifications (Note: num. classifications == num. random subsamples)
    - Keep predictions that do not contaminate previous bins:
9. Add filtered predictions to clusters
    - Re-select training_data with newly added contigs
10. Check base termination case:
    base case is when no filtered predictions are available
    if base case: return table
    else: jump back to 2.

"""


import logging
import random
import warnings

import numpy as np
import pandas as pd

from collections import namedtuple
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.tree import DecisionTreeClassifier

from autometa.common.markers import Markers


logger = logging.getLogger(__name__)

# Disable pandas and numpy warnings
pd.options.mode.chained_assignment = None
# https://stackoverflow.com/a/46721064/13118765
warnings.simplefilter(action="ignore", category=FutureWarning)

BinData = namedtuple("BinData", ["clustered", "unclustered"])
TrainingData = namedtuple("TrainingData", ["features", "target", "target_names"])


def get_taxa_features(fpath, dimensions=None):
    """Get a one-hot encoding of taxonomic information from `fpath`.

    Parameters
    ----------
    fpath : str
        Path to taxonomy table

    Returns
    -------
    pd.DataFrame
        One-hot encoded table of taxonomic features.
        index=contig, prefixes=['phylum','class','order','family','genus','species']
        columns will be prefixed with one of above and contain unique names for prefix-specific taxa.
    """
    df = pd.read_csv(fpath, sep="\t", index_col="contig")
    cols = ["phylum", "class", "order", "family", "genus", "species"]
    encoded_df = pd.get_dummies(df[cols], columns=cols, prefix=cols)
    # encoded_df = pd.get_dummies(df)
    if dimensions:
        pca = PCA(dimensions)
        X = encoded_df.to_numpy()
        X_fit = pca.fit_transform(X)
        encoded_df = pd.DataFrame(X_fit, index=encoded_df.index)
    return encoded_df


def get_kmer_features(fpath, dimensions=None):
    """Get k-mer features from normalized k-mer frequencies reduced by PCA to `dimensions`

    Parameters
    ----------
    fpath : str
        Path to normalized k-mer frequencies table
    dimensions : int, optional
        number of principal components to reduce k-mer frequencies down to, by default will not reduce by PCA

    Returns
    -------
    pd.DataFrame
        shape=(n_contigs, dimensions), index=contig, cols=range([0, dimensions])
    """
    df = pd.read_csv(fpath, sep="\t", index_col="contig")
    # logger.debug(f"PCA dimension reduction (dimensions={dimensions})")
    if dimensions:
        pca = PCA(n_components=dimensions)
        X = df.to_numpy()
        X_fit = pca.fit_transform(X)
        df = pd.DataFrame(X_fit, index=df.index)
    return df


def get_features(
    kmers,
    coverage,
    annotations=[],
    taxonomy=None,
    kmer_dimensions=50,
    taxa_dimensions=10,
):
    """Retrieve `dimensions` principal components from `kmers` then merge additional annotations to
    create a master dataframe of contig features.

    Parameters
    ----------
    kmers : str
        Path to normalized k-mer frequencies table which will be provided to PCA(dimensions)
    coverage: str
        Path to contig coverage calculations table
    annotations : list, optional
        Paths to additional annotations to add to features, by default []
    taxonomy : str, optional
        Path to taxonomic features, by default None
    dimensions : int, optional
        Number of principal components to retrieve from normalized k-mer frequencies, by default 50

    Returns
    -------
    pd.DataFrame
        index=contig, cols=[0-dimension, coverage, ...]
        additional columns would correspond to any features provided by `taxonomy` or `annotations`
    """
    df = get_kmer_features(fpath=kmers, dimensions=kmer_dimensions)
    # annotations is here in case you would like to add additional annotations as features
    annotations.append(coverage)
    for fpath in annotations:
        dff = pd.read_csv(fpath, sep="\t", index_col="contig")
        df = pd.merge(df, dff, how="inner", left_index=True, right_index=True)
    if taxonomy:
        taxa_df = get_taxa_features(taxonomy, taxa_dimensions)
        df = pd.merge(df, taxa_df, how="inner", left_index=True, right_index=True)

    omit_cols = {"cluster", "reference_genome", "completeness", "purity"}
    feature_cols = [col for col in df.columns if col not in omit_cols]
    return df[feature_cols]


def get_labels(df):
    """Retrieve clustering labels from `df`.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe containing either `reference_genome` or `cluster` annotations.

    Returns
    -------
    namedtuple("Labels", ['target':binary matrix of contig to label assignments, 'target_names': labels associated with binary matrix indices])
        `reference_genome` or `cluster` labels

    Raises
    ------
    ValueError
        Either `reference_genome` or `cluster` are not in `df.columns`
    """
    if "reference_genome" in df.columns:
        df = df[df["reference_genome"] != "misassembled"]
        target = pd.get_dummies(df["reference_genome"])
    elif "cluster" in df.columns:
        target = pd.get_dummies(df[df.cluster.notnull()]["cluster"])
    else:
        bin_cols = ["reference_genome", "cluster"]
        raise ValueError(f"{df.columns} does not contain one of bin_cols: {bin_cols}")
    Labels = namedtuple("labels", ["target", "target_names"])
    return Labels(target=target, target_names=target.columns.tolist())


def split_clustered_unclustered(bin_df):
    """Split `bin_df` into clustered and unclustered contigs

    Parameters
    ----------
    bin_df : pd.DataFrame
        dataframe containing 'cluster' col to split null and not null

    Returns
    -------
    namedtuple("BinData", ['clustered':pd.DataFrame, 'unclustered':pd.DataFrame])
        split dataframes of clustered and unclustered contigs
    """
    unclustered_df = bin_df.loc[bin_df.cluster.isnull()]
    clustered_df = bin_df.loc[bin_df.cluster.notnull()]
    return BinData(clustered=clustered_df, unclustered=unclustered_df)


def get_clustered_markers(bin_data, markers_df):
    """Retrieve clustered marker-containing contigs`

    Parameters
    ----------
    bin_data : namedtuple("BinData", ['clustered':pd.DataFrame,'unclustered':pd.DataFrame])
        binning data split between clustered/unclustered contigs
    markers_df : pd.DataFrame
        Dataframe retrieved from :func: Markers.load(fpath, format='wide')

    Returns
    -------
    namedtuple("BinData", ['clustered':pd.DataFrame,'unclustered':pd.DataFrame])
        binning data split between clustered/unclustered contigs where clustered is subset to
        only marker-containing contigs.
    """
    clustered_markers_index = bin_data.clustered.index.isin(markers_df.index)
    clustered_markers_df = bin_data.clustered.loc[clustered_markers_index]
    return BinData(clustered=clustered_markers_df, unclustered=bin_data.unclustered)


def split_and_subset_train_data(bin_data, features, labels):
    """Subset features and labels using split between clustered/unclustered contigs.

    Note
    ----
    namedtuple("TrainingData", ["features":pd.DataFrame, "target":pd.DataFrame, "target_names":list])

    Parameters
    ----------
    bin_data : namedtuple("BinData", ['clustered':pd.DataFrame,'unclustered':pd.DataFrame])
        Bin data split by clustered/unclustered contigs
    features : pd.DataFrame
        Contig features to be used for training a classifier and predicting classifications
    labels : pd.DataFrame
        Labels corresponding to classifier predictions

    Returns
    -------
    namedtuple("BinData", ['clustered':namedtuple("TrainingData"),'unclustered':namedtuple("TrainingData")]
        Features and bin labels split and subset by clustered/unclustered contigs
    """
    clustered_features_index = features.index.isin(bin_data.clustered.index)
    clustered_features = features.loc[clustered_features_index]
    clustered_labels_index = labels.target.index.isin(bin_data.clustered.index)
    clustered_labels = labels.target.loc[clustered_labels_index]
    clustered_data = TrainingData(
        features=clustered_features,
        target=clustered_labels,
        target_names=clustered_labels.columns.tolist(),
    )

    unclustered_features_index = features.index.isin(bin_data.unclustered.index)
    unclustered_features = features.loc[unclustered_features_index]
    unclustered_labels_index = labels.target.index.isin(bin_data.unclustered.index)
    unclustered_labels = labels.target.loc[unclustered_labels_index]
    unclustered_data = TrainingData(
        features=unclustered_features,
        target=unclustered_labels,
        target_names=unclustered_labels.columns.tolist(),
    )
    return BinData(clustered=clustered_data, unclustered=unclustered_data)


def get_decision_tree_predictions(
    X, y, y_test, bin_data, num_classifications=10, seed=42
):
    """Get predictions using DecisionTreeClassifier.

    Parameters
    ----------
    X : numpy.ndarray shape=(n_contigs, n_features)
        Clustered contig features to train the classifier.
    y : numpy.ndarray, shape=(n_labels, )
        Clustered contig labels to train the classifier.
    y_test : numpy.ndarray, shape=(n_contigs, n_features)
        Unclustered contig features for input to ``classifier.predict(y_test)``
    bin_data : namedtuple("BinData", ['clustered':namedtuple("TrainingData"),'unclustered':namedtuple("TrainingData")]
        Features and bin labels split and subset by clustered/unclustered contigs.
    num_classifications : int, optional
        number of classifications to perform, by default 10

    Returns
    -------
    pd.DataFrame
        Summed predictions by `num_classifications`. index=contig, columns=[label, label, ...]
    """
    predictions = []
    for i in range(num_classifications):
        X_train, _, y_train, _ = train_test_split(X, y, test_size=0.5, random_state=i)
        clf = DecisionTreeClassifier(random_state=np.random.RandomState(seed))
        clf.fit(X_train, y_train)
        prediction = clf.predict(y_test)
        predictions.append(prediction)

    summed_predictions = np.sum(predictions, axis=0)
    return pd.DataFrame(
        summed_predictions,
        index=bin_data.unclustered.features.index,
        columns=bin_data.clustered.target_names,
    )


def get_random_forest_predictions(X, y, y_test, bin_data, num_estimators=10, seed=42):
    """Retrieve predictions using RandomForestClassifier.

    Note
    ----
    By passing in `num_estimators` and `bootstrap=True` as parameters to this classifier,
    the mean probabilities across `num_estimators` predictions will be calculated and the highest
    class with the highest probability will be returned.

    Parameters
    ----------
    X : numpy.ndarray shape=(n_contigs, n_features)
        Clustered contig features to train the classifier.
    y : numpy.ndarray, shape=(n_labels, )
        Clustered contig labels to train the classifier.
    y_test : numpy.ndarray, shape=(n_contigs, n_features)
        Unclustered contig features for input to ``classifier.predict(y_test)``
    bin_data : namedtuple("BinData", ['clustered':namedtuple("TrainingData"),'unclustered':namedtuple("TrainingData")]
        Features and bin labels split and subset by clustered/unclustered contigs.
    num_estimators : int, optional
        number of estimators to construct for training and predictions, by default 10

    Returns
    -------
    pd.DataFrame
        Binary matrix of predictions from `num_estimators`. index=contig, columns=[label, label, ...]
    """
    clf = RandomForestClassifier(
        max_samples=0.5,
        n_estimators=num_estimators,
        bootstrap=True,
        random_state=np.random.RandomState(seed),
    )
    clf.fit(X, y)
    prediction = clf.predict(y_test)
    return pd.DataFrame(
        prediction,
        index=bin_data.unclustered.features.index,
        columns=bin_data.clustered.target_names,
    )


def get_confidence_filtered_predictions(
    bin_data,
    num_classifications=10,
    confidence=1.0,
    classifier="decision_tree",
    seed=42,
):
    """Filter classifier predictions by confidence threshold.

    Parameters
    ----------
    bin_data : namedtuple("BinData", ['clustered':namedtuple("TrainingData"),'unclustered':namedtuple("TrainingData")]
        Features and bin labels split and subset by clustered/unclustered contigs.
    num_classifications : num, optional
        Number of classifications to perform on each round of predictions to determine classifier confidence, by default 10
    confidence : float, optional
        Threshold to keep/discard classifier predictions (num. consistent classifications / num. classifications), by default 1.0
    classifier : str, optional
        classifier to use for predictions, by default "decision_tree"
        choices include "decision_tree" and "random_forest".

    Returns
    -------
    pd.DataFrame
        Confidence filtered predictions dataframe, index=contig, col='cluster'

    Raises
    ------
    NotImplementedError
        Provided `classifier` is not implemented.
    """
    X = bin_data.clustered.features.to_numpy()
    y = bin_data.clustered.target.to_numpy()
    y_test = bin_data.unclustered.features.to_numpy()

    logger.debug(
        f"getting predictions from {X.shape[0]} contigs with {X.shape[1]} features for {y_test.shape[0]} unclustered contigs"
    )

    if classifier == "decision_tree":
        df = get_decision_tree_predictions(
            X, y, y_test, bin_data, num_classifications=num_classifications, seed=seed
        )
    elif classifier == "random_forest":
        df = get_random_forest_predictions(
            X, y, y_test, bin_data, num_estimators=num_classifications, seed=seed
        )
    else:
        raise NotImplementedError(classifier)

    confidence_threshold = num_classifications * confidence
    df = df.loc[df.max(axis=1) >= confidence_threshold]
    filtered_predictions = df.idxmax(axis=1)
    filtered_predictions.name = "cluster"
    return filtered_predictions.to_frame()


def filter_contaminating_predictions(predictions_df, markers_df, bin_df):
    """Filter contigs that would cause contamination in predicted bin.

    Parameters
    ----------
    predictions_df : pd.DataFrame
        classifier predictions, index=contig, col='cluster'
    markers_df : pd.DataFrame
        markers retrieved from Markers.load(fpath, format='wide')
    bin_df : pd.DataFrame
        Binning dataframe with which to compare classifications.

    Returns
    -------
    pd.DataFrame
        contamination filtered predictions, index=contig, col='cluster'
    """
    for cluster, dff in bin_df.groupby("cluster"):
        prediction_index = predictions_df[predictions_df.cluster == cluster].index
        if prediction_index.empty:
            # No reason to perform calculations if no predictions exist for current cluster
            continue

        old_markers = markers_df.loc[markers_df.index.isin(dff.index)]
        old_marker_counts = old_markers.sum()

        prediction_markers_index = markers_df.loc[
            markers_df.index.isin(prediction_index)
        ].index
        new_markers_index = prediction_markers_index.join(
            old_markers.index, how="outer"
        )
        new_marker_counts = markers_df.loc[new_markers_index].sum()

        old_num_single_copy_markers = old_marker_counts[old_marker_counts == 1].count()
        new_num_single_copy_markers = new_marker_counts[new_marker_counts == 1].count()
        if new_num_single_copy_markers >= old_num_single_copy_markers:
            # We are not adding contamination, so keep the contig predictions
            continue

        # This means we would be introducing contamination so we drop the contigs from the predictions
        predictions_df.drop(index=prediction_index, inplace=True)

    return predictions_df


def add_predictions(bin_df, predictions_df):
    """Add classifier predictions to existing binning annotations.

    Parameters
    ----------
    bin_df : pd.DataFrame
        Existing binning dataframe to be updated
    predictions_df : pd.DataFrame
        classifier predictions to add into `bin_df`

    Returns
    -------
    pd.DataFrame
        Updated `bin_df`
    """
    for cluster, dff in bin_df.groupby("cluster"):
        prediction_index = predictions_df[predictions_df.cluster == cluster].index
        if prediction_index.empty:
            # No reason to perform calculations if no predictions exist for current cluster
            continue
        bin_df.loc[prediction_index, "cluster"] = cluster
    return bin_df


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Recruit unclustered contigs given metagenome annotations and Autometa binning results."
        " Note: All tables must contain a 'contig' column to be used as the unique table index)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("kmers", help="Path to normalized kmer frequencies table.")
    parser.add_argument("coverage", help="Path to coverage table.")
    parser.add_argument(
        "assignments",
        help="Path to contig to genome bin assignments."
        "(May be autometa binning output [col='cluster']"
        "or ground truth bin assignments [col='reference_genome'])",
    )
    parser.add_argument("markers", help="Path to domain-specific markers table.")
    parser.add_argument("out", help="Path to output unclustered recruitment table.")
    parser.add_argument("--taxonomy", help="Path to taxonomy table.")
    parser.add_argument(
        "--taxa-dimensions",
        help="Num of dimensions to reduce taxonomy encodings",
        type=int,
    )
    parser.add_argument(
        "--additional-features",
        help="Path to additional features with which to add to classifier training data.",
        nargs="*",
        default=[],
    )
    parser.add_argument(
        "--confidence",
        help="Percent confidence to allow classification (confidence = num. consistent predictions/num. classifications)",
        default=1.0,
        type=float,
    )
    parser.add_argument(
        "--num-classifications",
        help="Num classifications for predicting/validating contig cluster recruitment",
        default=10,
        type=int,
    )
    parser.add_argument(
        "--classifier",
        help="classifier to use for recruitment of contigs",
        default="decision_tree",
        choices=["decision_tree", "random_forest"],
    )
    parser.add_argument(
        "--kmer-dimensions",
        help="Num of dimensions to reduce normalized k-mer frequencies",
        default=50,
        type=int,
    )
    parser.add_argument(
        "--seed",
        help="Seed to use for RandomState when initializing classifiers.",
        default=42,
        type=int,
    )
    args = parser.parse_args()

    features = get_features(
        kmers=args.kmers,
        coverage=args.coverage,
        annotations=args.additional_features,
        taxonomy=args.taxonomy,
        kmer_dimensions=args.kmer_dimensions,
        taxa_dimensions=args.taxa_dimensions,
    )
    bin_df = pd.read_csv(
        args.assignments, sep="\t", index_col="contig", usecols=["contig", "cluster"]
    )
    prev_num_unclustered = bin_df[bin_df.cluster.isnull()].shape[0]
    labels = get_labels(bin_df)
    markers_df = Markers.load(fpath=args.markers, format="wide")

    logger.debug(
        f"classifier={args.classifier}, seed={args.seed}, n.estimators={args.num_classifications}, confidence={args.confidence*100}%"
    )

    n_runs = 0
    while True:
        n_runs += 1

        bin_data = split_clustered_unclustered(bin_df=bin_df)

        bin_data = get_clustered_markers(bin_data=bin_data, markers_df=markers_df)

        labels = get_labels(bin_df)
        bin_data = split_and_subset_train_data(bin_data, features, labels)

        # Perform cross-validation with n. iterations (num. estimators)
        predictions_df = get_confidence_filtered_predictions(
            bin_data=bin_data,
            num_classifications=args.num_classifications,
            confidence=args.confidence,
            classifier=args.classifier,
            seed=args.seed,
        )
        predictions_df = filter_contaminating_predictions(
            predictions_df=predictions_df, markers_df=markers_df, bin_df=bin_df
        )

        # n_clustered = predictions_df.shape[0]
        # recruited_to_bins = ",".join(predictions_df.cluster.unique())
        # logger.debug(f"{n_clustered} contigs clustered to {recruited_to_bins}")

        if predictions_df.empty:
            break

        bin_df = add_predictions(bin_df=bin_df, predictions_df=predictions_df)

    now_num_unclustered = bin_df[bin_df.cluster.isnull()].shape[0]

    n_recruited = prev_num_unclustered - now_num_unclustered
    logger.info(
        f"unclustered {prev_num_unclustered} -> {now_num_unclustered} (recruited {n_recruited} contigs) in {n_runs} runs"
    )
    bin_df.to_csv(args.out, sep="\t", index=True, header=True)


if __name__ == "__main__":
    main()
