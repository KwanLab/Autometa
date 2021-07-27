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
import numpy as np
import sys

input_path = sys.argv[1]


def summarize_f1_stats(f1_table_path):
    f1_table = pd.read_csv(f1_table_path, sep="\t")
    f1_dict = {}
    for count, row in f1_table.iterrows():
        if row["ref_genome"] != "misassembled" and str(row["ref_genome"]) != "nan":
            f1_dict[row["ref_genome"]] = row["F1"]
    return (
        f1_table_path,
        round(sum(f1_dict.values()), 1),
        round(np.median(list(f1_dict.values())), 1),
    )


print(summarize_f1_stats(input_path))
