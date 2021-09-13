#!/usr/bin/env python

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

# Docker wrapper for ML_recruitment.py
# Aims to be invisible to the user (use in exactly the same way as ML_recruitment.py)

import argparse
import os
import subprocess


def run_command(command_string, stdout_path=None):
    # Function that checks if a command ran properly. If it didn't, then print an error message then quit
    if stdout_path:
        f = open(stdout_path, "w")
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print("ML_recruitment_docker.py: Error, the command:")
        print(command_string)
        print(("failed, with exit code " + str(exit_code)))
        exit(1)


parser = argparse.ArgumentParser(
    description="Recruit unclustered (or non-marker)\
	sequences with Machine Learning classification using clustered sequence\
	as training data. Features to train with include sequence coverage,\
	composition, and homology. Confidence is calculated using jackknife\
	cross-validation by randomly subsetting the training data n number of times."
)
parser.add_argument(
    "-t",
    "--contig_tab",
    metavar="<contig.tab>",
    help="Path to master contig table which includes initial clusters",
    required=True,
)
parser.add_argument(
    "-c",
    "--cluster_column",
    metavar="<column header>",
    help="Name of column containing initial cluster information",
    default="cluster",
)
parser.add_argument(
    "-p",
    "--processors",
    metavar="<int>",
    help="Number of processors to use",
    type=int,
    default=1,
)
parser.add_argument(
    "-r",
    "--recursive",
    help="If specified, will run classification \
	iteratively and refine traning data after each iteration.",
    action="store_true",
)
parser.add_argument(
    "-C",
    "--Confidence_cutoff",
    metavar="<int>",
    help="Confidence cutoff value\
	to use to keep ML-based predictions.",
    type=int,
    default=100,
)
parser.add_argument(
    "-u",
    "--unclustered_name",
    metavar="<unclustered name>",
    help="Name of unclustered group \
	in cluster column",
    default="unclustered",
)
parser.add_argument(
    "-n",
    "--num_iterations",
    metavar="<int>",
    help="Number of iterations for \
	jackknife cross-validation.",
    type=int,
    default=10,
)
parser.add_argument(
    "-m",
    "--k_mer_matrix",
    metavar="<k-mer.tab>",
    help="Path to k-mer_matrix file.",
    default="k-mer_matrix",
)
parser.add_argument(
    "-o",
    "--out_table",
    metavar="<output.tab>",
    help="Path to create output table with new column\
	for ML-recruited sequences.",
    required=True,
)
parser.add_argument(
    "-k",
    "--kingdom",
    metavar="<archaea|bacteria>",
    help="Kingdom to consider (archaea|bacteria)",
    choices=["bacteria", "archaea"],
    default="bacteria",
)
args = vars(parser.parse_args())

contig_tab_path = args["contig_tab"]
cluster_column = args["cluster_column"]
processors = args["processors"]
recursive = args["recursive"]
confidence_cutoff = args["Confidence_cutoff"]
unclustered_name = args["unclustered_name"]
num_iterations = args["num_iterations"]
k_mer_matrix_path = args["k_mer_matrix"]
out_table_path = args["out_table"]
kingdom = args["kingdom"]

# Do check to see if the contig table and the k_mer_matrix file exists
if not os.path.isfile(args["contig_tab"]):
    print(
        (
            "Error! Could not find contig table at the following path: "
            + args["contig_tab"]
        )
    )
    exit(1)

if not os.path.isfile(args["k_mer_matrix"]):
    print(
        (
            "Error! Could not find k-mer matrix file at the following path: "
            + args["k_mer_matrix"]
        )
    )
    exit(1)

# We will use the directory that will hold the output table as the output directory
out_table_path_absolute = os.path.abspath(out_table_path)
output_dir = "/".join(out_table_path_absolute.split("/")[:-1])
out_table_filename = out_table_path_absolute.split("/")[-1]

# We have to create the output dir if it doesn't exist
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

# If the contig table is not already in the output dir, we need to copy it there so that the docker container can see it
contig_tab_path_absolute = os.path.abspath(contig_tab_path)
contig_tab_directory = "/".join(contig_tab_path_absolute.split("/")[:-1])
contig_tab_filename = contig_tab_path_absolute.split("/")[-1]

if output_dir != contig_tab_directory:
    # This means we need to copy the contig table to the output directory
    run_command("cp " + contig_tab_path_absolute + " " + output_dir + "/")

# If k-mer matrix file is not already in the output dir, we need to copy it there so that the docker container can see it
k_mer_matrix_path_absolute = os.path.abspath(k_mer_matrix_path)
k_mer_matrix_directory = "/".join(k_mer_matrix_path_absolute.split("/")[:-1])
k_mer_matrix_filename = k_mer_matrix_path_absolute.split("/")[-1]

if output_dir != k_mer_matrix_directory:
    run_command("cp " + k_mer_matrix_path_absolute + " " + output_dir + "/")

# Construct ML_recruitment.py command to pass to the docker container
ML_recruitment_command = "ML_recruitment.py --contig_tab /output/{} --cluster_column {} --processors {} --Confidence_cutoff {} --unclustered_name {} --num_iterations {} --k_mer_matrix /output/{} --out_table /output/{} --kingdom {}".format(
    contig_tab_filename,
    cluster_column,
    processors,
    confidence_cutoff,
    unclustered_name,
    num_iterations,
    k_mer_matrix_filename,
    out_table_filename,
    kingdom,
)

if recursive:
    ML_recruitment_command = ML_recruitment_command + " --recursive"

# Construct Docker run command
docker_command = "docker run --volume {}:/output:rw --detach=false --rm jasonkwan/autometa:latest".format(
    output_dir
)

# Incorporate autometa command
docker_command = docker_command + " " + ML_recruitment_command

# Run docker container
run_command(docker_command)
