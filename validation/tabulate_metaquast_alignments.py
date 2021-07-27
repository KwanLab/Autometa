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

import subprocess
import argparse

parser = argparse.ArgumentParser(
    description="Script to tabulate reference genomes\
    from metaquast NUCMER alignments."
)
parser.add_argument("-p", "--asm_prefix", help="Assembly prefix name", required=True)
parser.add_argument(
    "-o",
    "--out",
    help="Tab delimited output table with reference genomes",
    default="metaquast_assigned_reference_genome.tab",
)
args = vars(parser.parse_args())

alignment_dict = {}
line_list = []
key_indices = []

asm_prefix = args["asm_prefix"]
all_alignments_file = "all_alignments_{}.tsv".format(asm_prefix)
alignments_file = "alignments_{}.tsv".format(asm_prefix)

with open(all_alignments_file) as infile:
    # line_list = [line for line in infile]
    # skip first line
    for count, line in enumerate(infile):
        if count > 0:
            line_list.append(line)
            adjusted_line_index = count - 1
            if "CONTIG" in line:
                key_indices.append(adjusted_line_index)

temp_list = []
# Skip header
for count, line in enumerate(line_list):
    if count in key_indices:
        split_line = line.rstrip().split("\t")
        contig = split_line[1]
        contig_status = split_line[3]
        # Add all of the lines to dictionary by CONTIG key (so that all info is captured)
        alignment_dict[contig] = [line for line in temp_list]
        # Then flush list for next round
        temp_list = []
    else:
        temp_list.append(line)

output_dict = {}
num_unaligned_contigs = 0
with open(args["out"], "w") as outfile:
    outfile.write("contig\treference_genome\n")
    for contig, info_list in list(alignment_dict.items()):
        output_dict[contig] = {}
        output_dict[contig]["misassembly_type"] = "NA"
        print(contig)
        contig_length = int(contig.split("_")[3])
        # If there's not detailed alignment info, check status in first file
        if len(info_list) == 0:
            aln_status = subprocess.check_output(
                "grep {} {}".format(contig, all_alignments_file), shell=True
            )
            if aln_status.split()[3] == "correct":
                # Alignment info (to ref genome)  stored in another file..
                reference_genome = (
                    subprocess.check_output(
                        "grep {} {}".format(contig, alignments_file), shell=True
                    )
                ).split()[0]
                # print reference_genome
                output_dict[contig]["reference_genome"] = reference_genome
                output_dict[contig]["id_to_reference_genome"] = 100.0
                output_dict[contig]["aln_length"] = contig_length
            else:
                print(("{} is unaligned!".format(contig)))
                output_dict[contig]["reference_genome"] = "NA"
                output_dict[contig]["id_to_reference_genome"] = "NA"
                output_dict[contig]["aln_length"] = "NA"
                num_unaligned_contigs += 1
                # continue

        else:
            for line in info_list:
                # Information on the ultimate mis-assembly type is alone one line and only appears once
                if len(line.split("\t")) == 1:
                    misassembly_type = line.rstrip()
                    if misassembly_type == "interspecies translocation":
                        output_dict[contig]["reference_genome"] = "misassembled"
                        # continue
                # Extra alignment info can be on multiple lines and always has eight values
                elif len(line.split("\t")) == 8:
                    # print line.rstrip()
                    line = line.split("\t")
                    reference_genome = line[4]
                    id_to_reference_genome = line[6]
                    aln_length = abs(int(line[3]) - int(line[2])) + 1
                    assert aln_length <= contig_length
                    if "reference_genome" not in output_dict[contig]:
                        output_dict[contig]["reference_genome"] = reference_genome
                        output_dict[contig][
                            "id_to_reference_genome"
                        ] = id_to_reference_genome
                        output_dict[contig]["aln_length"] = aln_length
                    # Set reference geneome to the alignment with the highest ID - .. might want to take len into consideration...
                    elif (
                        id_to_reference_genome
                        > output_dict[contig]["id_to_reference_genome"]
                    ):
                        output_dict[contig]["reference_genome"] = reference_genome
                        output_dict[contig][
                            "id_to_reference_genome"
                        ] = id_to_reference_genome
                        output_dict[contig]["aln_length"] = aln_length

        # print("--> {} {} {}".format(output_dict[contig]['reference_genome'].split(".")[0],output_dict[contig]['id_to_reference_genome'],output_dict[contig]['aln_length']))
        outfile.write(
            contig + "\t" + output_dict[contig]["reference_genome"].split(".")[0] + "\n"
        )

# Number of contigs in output table may be different from number of contigs in asm file.
# Based on unaligned_report.txt, seems that partially unaligned contigs with misassembly
# are absent from all_alignments_*.tsv file.
