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

import glob
import argparse
from os import path

parser = argparse.ArgumentParser(
    description="Script to build contig-cluster table from fasta files."
)
parser.add_argument(
    "-i", "--input_path", help="Path to directory with fasta files.", required=True
)
parser.add_argument(
    "-c",
    "--cluster_name",
    help="Name of cluster column (Default: cluster)",
    default="cluster",
)
parser.add_argument(
    "-x", "--extension", help="Fasta file extension (Default: fasta)", default="fasta"
)
parser.add_argument(
    "-o", "--outfile", help="Output table path", default="contig_bins.tab"
)
args = vars(parser.parse_args())

dir_path = path.realpath(args["input_path"])
fasta_regex = dir_path + "/*." + args["extension"]

with open(args["outfile"], "w") as outfile:
    outfile.write("contig\t{}\n".format(args["cluster_name"]))
    for fasta in glob.glob(fasta_regex):
        try:
            fasta_handle = open(fasta)
            for line in fasta_handle:
                if ">" in line:
                    contig_name = line.rstrip().lstrip(">")
                    bin_name = ".".join(path.basename(fasta).split(".")[:-1])
                    outfile.write(contig_name + "\t" + bin_name + "\n")
        except IOError as e:
            print(e)
            print(("Issue opening {}\nExiting...".format(fasta)))
