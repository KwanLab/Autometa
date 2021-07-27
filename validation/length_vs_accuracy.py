#!/usr/bin/env python2.7

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import argparse
import math

parser = argparse.ArgumentParser(
    description="Script to tabulate reference training table (training based on marker contigs)."
)
parser.add_argument(
    "-t", "--contig_tab", help="Name of master contig table file", required=True
)
# parser.add_argument('-o','--out_tab', help='Name of output table file', required=True)
args = vars(parser.parse_args())

# master_table = pd.read_csv("18-Aug-17_MIX_51_w_ML_recruitment_C0n10p2.tab", sep = '\t')
master_table = pd.read_csv(args["contig_tab"], sep="\t")

prediction_accuracy_dict = {}
unclustered_length = 0
correctly_classified_length = 0
training_data_length = 0
total_length_over10k = 0
for count, row in master_table.iterrows():
    length = row["length"]
    if length >= 10000 and row["reference_training"] == "unclustered":
        unclustered_length += length
        ML_prediction = row["ML_expanded_clustering"]
        if ML_prediction == row["reference_genome"]:
            accurate_prediction = True
            correctly_classified_length += length
        else:
            accurate_prediction = False
        length_cutoff = math.ceil(length / 10000.0) * 10000.0
        if not length_cutoff in prediction_accuracy_dict:
            prediction_accuracy_dict[length_cutoff] = [accurate_prediction]
        else:
            prediction_accuracy_dict[length_cutoff].append(accurate_prediction)
    elif length >= 10000 and row["reference_training"] != "unclustered":
        training_data_length += length
    if length >= 10000:
        total_length_over10k += length

total_accurate_predictions = 0
total_predictions = 0
print("len_cutoff\tpercent_accurate\tnum_predictions\tnum_accurate_predictions")
for length_cutoff in range(10000, 500000, 10000):
    number_of_predictions = 0
    rounded_length_cutoff = str(length_cutoff)[:-3]
    try:
        number_of_accurate_predictions = prediction_accuracy_dict[length_cutoff].count(
            True
        )
        total_accurate_predictions += number_of_accurate_predictions
        number_of_predictions = len(prediction_accuracy_dict[length_cutoff])
        total_predictions += number_of_predictions
        percent_accuracy = str(
            round(
                number_of_accurate_predictions / float(number_of_predictions) * 100, 2
            )
        )
        print(
            (
                rounded_length_cutoff
                + "\t"
                + percent_accuracy
                + "\t"
                + str(number_of_predictions)
                + "\t"
                + str(number_of_accurate_predictions)
            )
        )
    except KeyError:
        print((rounded_length_cutoff + "\t" + "NA" + "\t" + str(0) + "\t" "NA" + "\t"))

print(
    (
        "Initially unclustered length is {} and ML classified length is {}".format(
            unclustered_length, correctly_classified_length
        )
    )
)
test_data_percent = round(unclustered_length / float(total_length_over10k) * 100, 2)
training_data_percent = round(
    training_data_length / float(total_length_over10k) * 100, 2
)
average_prediction_accuracy = round(
    total_accurate_predictions / float(total_predictions) * 100, 2
)
# average_length_weighted_prediction_accuracy =
print(
    (
        "Training and test data represent {}% and {}% of the total contigs, respectively. Average prediction accuracy: {}".format(
            training_data_percent, test_data_percent, average_prediction_accuracy
        )
    )
)
