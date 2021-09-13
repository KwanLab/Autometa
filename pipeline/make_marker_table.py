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

import pandas as pd
import argparse
import subprocess
import os

# argument parser
parser = argparse.ArgumentParser(
    description='Script tabulate single copy markers \
	from a metagenome assembly. Dependencies: prodigal v2.6.2 (from "GoogleImport" branch), hhmscan (hmmer 3.1b2)'
)
parser.add_argument("-a", "--assembly", help="Input assembly file", required=True)
parser.add_argument(
    "-p", "--processors", help="Number of processors to use for hmmscan", default=1
)
parser.add_argument(
    "-c",
    "--cutoffs",
    help="Bacterial single copy hmm cutoffs as defined by Rinke et al. Default path is home directory.",
    default="~/Bacteria_single_copy_cutoffs.txt",
)
parser.add_argument(
    "-m",
    "--hmm",
    help="Bacteria_single_copy_cutoffs.hmm. Default path is home directory.",
    default="~/Bacteria_single_copy.hmm",
)
parser.add_argument(
    "-o",
    "--out",
    help="outfile.tab, three column table with contig, single copy PFAMS, and # of markers",
    required=False,
)
args = vars(parser.parse_args())

assembly = os.path.abspath(args["assembly"])


def get_contig_list(path_to_assembly):
    # Get list of spades contigs
    assembly_handle = open(path_to_assembly, "rU")
    contig_name_list = []
    for line in assembly_handle:
        if ">" in line:
            contig_name = line.rstrip("\n").split()[0][1:]
            contig_name_list.append(contig_name)
    assembly_handle.close()
    return contig_name_list


def run_prodigal(path_to_assembly):
    assembly_filename = path_to_assembly.split("/")[-1]
    # When "shell = True", need to give one string, not a list
    subprocess.call(
        " ".join(
            [
                "prodigal ",
                "-i " + path_to_assembly,
                "-a " + output_dir + "/" + assembly_filename + ".orfs.faa",
                "-p meta",
                "-m",
                "-o " + output_dir + "/" + assembly_filename + ".txt",
            ]
        ),
        shell=True,
    )
    return output_dir + "/" + assembly_filename + ".orfs.faa"


def run_hhmscan(path_to_prodigal_output, hmmdb):
    subprocess.call(
        "hmmscan --cpu {} --tblout {} {} {}".format(
            args["processors"],
            path_to_prodigal_output + ".hmm.tbl",
            hmmdb,
            path_to_prodigal_output,
        ),
        shell=True,
    )
    return path_to_prodigal_output + ".hmm.tbl"


output_dir = "/".join(os.path.abspath(args["out"]).split("/")[:-1])

prodigal_output = run_prodigal(assembly)
hmm_table_path = run_hhmscan(prodigal_output, args["hmm"])

hmm_table = pd.read_csv(
    hmm_table_path,
    sep="\s+",
    usecols=[1, 2, 5],
    skiprows=3,
    header=None,
    index_col=False,
    engine="python",
)
cutoffs_table = pd.read_csv(args["cutoffs"], sep="\s", engine="python", header=None)

# Search for contigs/ORFs that contain single copy PFAM domains that pass cutoff
# for loop to search for PFAM domains in hmm table column 1:
contig_ORFs_that_pass_cutoffs = {}
for index, PFAM_cutoffs_id in enumerate(cutoffs_table[0]):
    for count, PFAM_hmm_scan_id in enumerate(hmm_table[1]):
        # If the PFAMs domain match (cutoffs ID names format are PFAMXXXXX, not PFAMXXXXX.1)
        if str(PFAM_cutoffs_id) in str(PFAM_hmm_scan_id):
            # and hmm score value > cutoff value,
            if float(hmm_table[5][count]) > float(cutoffs_table[1][index]):
                # make dictionary of lists for PFAMs, where PFAM is key and contigs populate the list
                if PFAM_hmm_scan_id not in contig_ORFs_that_pass_cutoffs:
                    contig_ORFs_that_pass_cutoffs[PFAM_hmm_scan_id] = []
                    contig_ORFs_that_pass_cutoffs[PFAM_hmm_scan_id].append(
                        hmm_table[2][count]
                    )
                else:
                    contig_ORFs_that_pass_cutoffs[PFAM_hmm_scan_id].append(
                        hmm_table[2][count]
                    )

# Make a dictionary of dictionaries with contigs that have single copy genes (list PFAMS, for final table count length of list),
# contig length, contig GC, contig len, passecd PFAM domains
# write out tab-delimited table

if args["out"] != None:
    outfile_handle = args["out"]
else:
    outfile_handle = assembly + ".marker.tab"
with open(outfile_handle, "w") as outfile:
    outfile.write(
        "contig" + "\t" + "single_copy_PFAMs" + "\t" + "num_single_copies" + "\n"
    )
    contig_dictionary = {}
    for count, contig in enumerate(get_contig_list(assembly)):
        contig_dictionary[contig] = {}
        contig_dictionary[contig]["single_copy_PFAMs"] = []
        contig_dictionary[contig]["num_single_copies"] = 0
        for PFAM_key, contigs in list(contig_ORFs_that_pass_cutoffs.items()):
            for item in contigs:
                if str(contig) in item:
                    contig_dictionary[contig]["single_copy_PFAMs"].append(PFAM_key)
        contig_dictionary[contig]["num_single_copies"] = len(
            contig_dictionary[contig]["single_copy_PFAMs"]
        )

        if len(contig_dictionary[contig]["single_copy_PFAMs"]) > 0:
            outfile.write(
                str(contig)
                + "\t"
                + str(",".join(contig_dictionary[contig]["single_copy_PFAMs"]))
                + "\t"
                + str(contig_dictionary[contig]["num_single_copies"])
                + "\n"
            )
        else:
            outfile.write(
                str(contig)
                + "\t"
                + "NA"
                + "\t"
                + str(contig_dictionary[contig]["num_single_copies"])
                + "\n"
            )

print("\nDone!")
