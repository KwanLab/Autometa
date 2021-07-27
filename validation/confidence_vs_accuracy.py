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
import math
import argparse

parser = argparse.ArgumentParser(
    description="Script to tabulate confidence vs accuracy."
)
parser.add_argument(
    "-t", "--contig_tab", help="Name of master contig table file", required=True
)
parser.add_argument(
    "-s",
    "--ML_recruitment_stdout",
    help="Name of master contig table file",
    required=True,
)
args = vars(parser.parse_args())

reference_genome_dict = {}
master_table = pd.read_csv(args["contig_tab"], sep="\t")
for count, row in master_table.iterrows():
    contig = row["contig"]
    reference_genome = row["reference_genome"]
    reference_genome_dict[contig] = reference_genome

confidence_dict = {}
bool_confidence_dict = {}
contig_dict = {}
# Ex. input "18-Aug-17_MIX_51_master_table_with_ML_recruitment.stdout"
with open(args["ML_recruitment_stdout"]) as infile:
    for line in infile:
        if "ML predictions and jackknife confidence for contig " in line:
            contig = line.split()[7][:-1]
            length = int(contig.split("_")[3])
            prediction = line.split()[8].split(",")[0]
            confidence = line.split()[8].split(",")[1]
            if prediction == reference_genome_dict[contig]:
                length_weighted_accuracy = length
                accurate_prediction = True
            elif prediction != reference_genome_dict[contig]:
                length_weighted_accuracy = -length
                accurate_prediction = False
            if not confidence in confidence_dict:
                confidence_dict[confidence] = [length_weighted_accuracy]
                bool_confidence_dict[confidence] = [accurate_prediction]
            else:
                confidence_dict[confidence].append(length_weighted_accuracy)
                bool_confidence_dict[confidence].append(accurate_prediction)
print(
    "confidence_interval\tsum_length_weighted_accuracy\tnum_accurate_predictions\tnum_inaccurate_predictions"
)
for confidence, length_weighted_accuracy_list in list(confidence_dict.items()):
    # sum_length_in_Mbp = math.ceil(round(sum(length_weighted_accuracy_list),4) / 1000000.0000) * 1000000.0000
    sum_length_in_Mbp = round(sum(length_weighted_accuracy_list) / 1000000, 4)
    # rounded_length_sum = str(sum_length_in_10kbp)[:-3]
    num_accurate_predictions = str(bool_confidence_dict[confidence].count(True))
    num_inaccurate_predictions = str(-(bool_confidence_dict[confidence].count(False)))
    print(
        (
            str(confidence).split(".")[0]
            + "\t"
            + str(sum_length_in_Mbp)
            + "\t"
            + num_accurate_predictions
            + "\t"
            + num_inaccurate_predictions
        )
    )
