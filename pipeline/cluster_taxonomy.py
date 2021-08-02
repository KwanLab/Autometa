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

# Program that determines the taxonomy of clusters called by dbcan
# Uses contig taxonomy information in the same way that add_contig_taxonomy.py uses protein taxonomy
# (except this program uses contig length as weighting rather than 1 protein 1 vote)
# Algorithm:
# In descending order of species, genus, family, order, class, phylum, superkingdom:
#    In descending order of votes:
#        If classification shares common ancestry with majority of other proteins, accept result
#    If no result, move up to next taxonomic level

import sys
from time import *
from tqdm import *
import subprocess
import argparse
import pprint

pp = pprint.PrettyPrinter(indent=4)

rank_priority = [
    "species",
    "genus",
    "family",
    "order",
    "class",
    "phylum",
    "superkingdom",
    "root",
]
canonical_ranks = {
    "superkingdom": 1,
    "phylum": 1,
    "class": 1,
    "order": 1,
    "family": 1,
    "genus": 1,
    "species": 1,
}


def isConsistentWithOtherOrfs(taxid, rank, contigDictionary, taxidDictionary):
    # Function that determines for a given taxid, whether the majority of proteins
    # in a contig, with rank equal to or above the given rank, are common ancestors of
    # the taxid.  If the majority are, this function returns True, otherwise it returns
    # False

    # First we make a modified rank_priority list that only includes the current rank and above
    ranks_to_consider = None
    for i, rankName in enumerate(rank_priority):
        if rankName == rank:
            ranks_to_consider = rank_priority[i:]
            break

    # Now we total up the consistent and inconsistent ORFs

    consistentTotal = 0
    inconsistentTotal = 0

    for rankName in ranks_to_consider:
        if rankName in contigDictionary:
            for current_taxid in contigDictionary[rankName]:
                if isCommonAncestor(current_taxid, taxid, taxidDictionary):
                    consistentTotal += contigDictionary[rankName][current_taxid]
                else:
                    inconsistentTotal += contigDictionary[rankName][current_taxid]

    if consistentTotal > inconsistentTotal:
        return True
    else:
        return False


def isCommonAncestor(potentialParentTaxid, childTaxid, taxidDictionary):
    current_taxid = childTaxid
    while int(current_taxid) != 1:
        if potentialParentTaxid == current_taxid:
            return True
        current_taxid = taxidDictionary[current_taxid]["parent"]
    return False


def lowest_majority(contigDictionary, taxidDictionary):
    taxid_totals = (
        {}
    )  # Dictionary of dictionary, keyed by rank then by taxid, holds totals accounting for whole taxid paths

    for rank in rank_priority:
        if rank in contigDictionary:
            ranks_to_consider = None
            for i, rankName in enumerate(rank_priority):
                if rankName == rank:
                    ranks_to_consider = rank_priority[i:]
                    break

            for taxid in contigDictionary[rank]:
                # Make a dictionary to total the number of canonical ranks hit while traversing the path
                # - so that we can add 'unclassified' to any that don't exist
                # Later we need to make sure that 'unclassified' doesn't ever win
                ranks_in_path = {}
                for rank_to_consider in ranks_to_consider:
                    ranks_in_path[rank_to_consider] = 0

                # We need to add to taxid_totals for each taxid in the tax_path
                current_taxid = taxid
                current_rank = rank
                while int(current_taxid) != 1:
                    if current_rank in canonical_ranks:
                        ranks_in_path[current_rank] += 1
                        if current_rank in taxid_totals:
                            if current_taxid in taxid_totals[current_rank]:
                                taxid_totals[current_rank][current_taxid] += 1
                            else:
                                taxid_totals[current_rank][current_taxid] = 1
                        else:
                            taxid_totals[current_rank] = {current_taxid: 1}
                    current_taxid = taxidDictionary[current_taxid]["parent"]
                    current_rank = taxidDictionary[current_taxid]["rank"]

                # Now go through ranks_in_path. Where total = 0, add 'unclassified'
                for rank_to_consider in ranks_to_consider:
                    if ranks_in_path[rank_to_consider] == 0:
                        if rank_to_consider in taxid_totals:
                            if "unclassified" in taxid_totals[rank_to_consider]:
                                taxid_totals[rank_to_consider]["unclassified"] += 1
                            else:
                                taxid_totals[rank_to_consider]["unclassified"] = 1
                        else:
                            taxid_totals[rank_to_consider] = {"unclassified": 1}

    # If there are any gaps in the taxonomy paths for any of the proteins in the contig,
    # we need to add 'unclassified' to the relevant canonical taxonomic rank.
    # However, we must never allow 'unclassified' to win! (That just won't really tell us anything)

    # Now we need to determine which is the first level to have a majority
    for rank in rank_priority:
        total_votes = 0
        taxid_leader = None
        taxid_leader_votes = 0
        if rank in taxid_totals:
            for taxid in taxid_totals[rank]:
                total_votes += taxid_totals[rank][taxid]
                if taxid_totals[rank][taxid] > taxid_leader_votes:
                    taxid_leader = taxid
                    taxid_leader_votes = taxid_totals[rank][taxid]
            majority_threshold = float(total_votes) / 2
            if (
                taxid_leader_votes > majority_threshold
                and taxid_leader != "unclassified"
            ):
                return taxid_leader

    # Just in case
    return 1


parser = argparse.ArgumentParser(
    description="Summarize the taxonomy of clusters in a table"
)
parser.add_argument("-t", "--contig_tab", help="Master contig table", required=True)
parser.add_argument(
    "-c", "--cluster_column", help="Name of column for cluster", default="cluster"
)
parser.add_argument(
    "-x",
    "--taxdump_path",
    help="Path to taxdump directory (downloaded from NCBI taxonomy)",
    required=True,
)
parser.add_argument("-o", "--out_table", help="Output table path", required=True)
parser.add_argument(
    "-s", "--single_genome", help="Specifies single genome mode", action="store_true"
)

args = vars(parser.parse_args())

contig_table_path = args["contig_tab"]
taxdump_dir_path = args["taxdump_path"]
output_file_path = args["out_table"]
cluster_column_name = args["cluster_column"]
single_genome_mode = args["single_genome"]

# Process NCBI taxdump files
names_dmp_path = taxdump_dir_path + "/names.dmp"
nodes_dmp_path = taxdump_dir_path + "/nodes.dmp"

taxids = {}
print(strftime("%Y-%m-%d %H:%M:%S") + " Processing taxid names")
wc_output = subprocess.check_output(["wc", "-l", names_dmp_path], text=True)
wc_list = wc_output.split()
number_of_lines = int(wc_list[0])

with open(names_dmp_path) as names_dmp:
    for line in tqdm(names_dmp, total=number_of_lines):
        line_list = line.rstrip("\n").split("|")
        # Remove trailing and leading spaces
        for i, value in enumerate(line_list):
            line_list[i] = value.strip()

        # Only add scientific name entries
        if line_list[3] == "scientific name":
            # line_list[1] = line_list[1].replace(' ', '_')
            taxids[int(float(line_list[0]))] = {"name": line_list[1]}

print(strftime("%Y-%m-%d %H:%M:%S") + " Processing taxid nodes")
wc_output = subprocess.check_output(["wc", "-l", nodes_dmp_path], text=True)
wc_list = wc_output.split()
number_of_lines = int(wc_list[0])

with open(nodes_dmp_path) as nodes_dmp:
    for line in tqdm(nodes_dmp, total=number_of_lines):
        line_list = line.rstrip("\n").split("|")
        # Remove trailing and leading spaces
        for i, value in enumerate(line_list):
            line_list[i] = value.strip()

        taxids[int(float(line_list[0]))]["parent"] = int(float(line_list[1]))
        taxids[int(float(line_list[0]))]["rank"] = line_list[2]

print(strftime("%Y-%m-%d %H:%M:%S") + " Parsing taxonomy table")
contig_classifications = {}
length_of_contigs = {}

# Work out number of lines in file
wc_output = subprocess.check_output(["wc", "-l", contig_table_path], text=True)
wc_list = wc_output.split()
number_of_lines = int(wc_list[0])

# Determine contig, length and cluster indexes
contig_table = open(contig_table_path, "r")
contig_table_lines = contig_table.read().splitlines()
contig_table.close
contig_table_first_line_list = contig_table_lines[0].split("\t")
column_count = {}
contig_index = None
cluster_index = None
length_index = None
taxid_index = None
for i, heading in enumerate(contig_table_first_line_list):
    if heading == "contig":
        contig_index = i
    if heading == cluster_column_name:
        cluster_index = i
    if heading == "length":
        length_index = i
    if heading == "taxid":
        taxid_index = i

    if heading in column_count:
        column_count[heading] += 1
    else:
        column_count[heading] = 1

if contig_index is None:
    print('Error, could not find a "contig" column in ' + contig_table_path)
    sys.exit(2)
if cluster_index is None and not single_genome_mode:
    print('Error, could not find a "cluster" column in ' + contig_table_path)
    sys.exit(2)
if length_index is None:
    print('Error, could not find a "length" column in  ' + contig_table_path)
    sys.exit(2)
if taxid_index is None:
    print('Error, could not find a "taxid" column in ' + contig_table_path)
    sys.exit(2)
if column_count["contig"] > 1:
    print('Error, there is more than one "contig" column in ' + contig_table_path)
    sys.exit(2)
if not single_genome_mode and column_count["cluster"] > 1:
    print('Error, there is more than one "cluster" column in ' + contig_table_path)
    sys.exit(2)
if column_count["length"] > 1:
    print('Error, there is more than one "length" column in ' + contig_table_path)
    sys.exit(2)
if column_count["taxid"] > 1:
    print('Error, there is more than one "taxid" column in ' + contig_table_path)

for i, line in enumerate(tqdm(contig_table_lines, total=number_of_lines)):
    if i > 0:
        line_list = line.rstrip("\n").split("\t")
        contig = line_list[contig_index]
        if single_genome_mode:
            cluster = "unclustered"
        else:
            cluster = line_list[cluster_index]
        length = line_list[length_index]
        try:
            taxid = int(float(line_list[taxid_index]))
        except ValueError:
            taxid = "unclassified"

        taxRank = None
        if taxid == "unclassified":
            taxid = 1
            taxRank = "root"
        elif taxid not in taxids:
            # This happens sometimes when the taxid database and the NR database are not in sync
            taxid = 1
            taxRank = "root"
        else:
            taxRank = taxids[taxid]["rank"]

        # Now get the taxid of the next canonical rank (if applicable)
        if taxRank == "no rank":
            taxid = 1
            taxRank = "root"

        if taxid != 1:
            while taxRank not in rank_priority:
                taxid = taxids[taxid]["parent"]
                taxRank = taxids[taxid]["rank"]

        # Keep running total of taxids for each cluster
        if cluster not in contig_classifications:
            contig_classifications[cluster] = {}
        if taxRank not in contig_classifications[cluster]:
            contig_classifications[cluster][taxRank] = {}

        if taxid not in contig_classifications[cluster][taxRank]:
            contig_classifications[cluster][taxRank][taxid] = int(
                float(length)
            )  # We need this because directly inting a string in the form '16079.0' doesn't work
        else:
            contig_classifications[cluster][taxRank][taxid] += int(float(length))

        # Count length of contigs per cluster
        if contig in length_of_contigs:
            length_of_contigs[contig] += int(float(length))
        else:
            length_of_contigs[contig] = int(float(length))

print(strftime("%Y-%m-%d %H:%M:%S") + " Ranking taxids")
top_taxids = {}
total_clusters = len(contig_classifications)

for cluster in tqdm(contig_classifications, total=total_clusters):
    acceptedTaxid = None
    for rank in rank_priority:
        if acceptedTaxid is not None:
            break
        # Order in descending order of votes
        if rank in contig_classifications[cluster]:
            ordered_taxids = sorted(
                contig_classifications[cluster][rank],
                key=contig_classifications[cluster][rank].__getitem__,
                reverse=True,
            )
            # sys.exit()
            for taxid in ordered_taxids:
                if isConsistentWithOtherOrfs(
                    taxid, rank, contig_classifications[cluster], taxids
                ):
                    acceptedTaxid = taxid
                    break

    # If acceptedTaxid is still None at this point, there was some kind of draw, so we need to find the lowest taxonomic level where there is a
    # majority
    if acceptedTaxid is None:
        acceptedTaxid = lowest_majority(contig_classifications[cluster], taxids)

    top_taxids[cluster] = acceptedTaxid

print(strftime("%Y-%m-%d %H:%M:%S") + " Resolving taxon paths")
taxon_paths = {}
# Dictionary of dictionaries, keyed by contig then rank, contains the taxon names
for cluster in tqdm(top_taxids, total=len(top_taxids)):
    taxon_paths[cluster] = {}
    current_taxid = top_taxids[cluster]

    while int(current_taxid) != 1:
        current_rank = taxids[current_taxid]["rank"]
        if current_rank in canonical_ranks:
            taxon_paths[cluster][current_rank] = taxids[current_taxid]["name"]
        current_taxid = taxids[current_taxid]["parent"]

    for rank in rank_priority:
        if rank not in taxon_paths[cluster]:
            taxon_paths[cluster][rank] = "unclassified"

print(strftime("%Y-%m-%d %H:%M:%S") + " Writing table")
output_table = open(output_file_path, "w")
if single_genome_mode:
    output_table.write("kingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\ttaxid\n")
else:
    output_table.write(
        "cluster\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\ttaxid\n"
    )

for cluster in taxon_paths:
    if single_genome_mode:
        output_table.write(
            str(taxon_paths[cluster]["superkingdom"])
            + "\t"
            + str(taxon_paths[cluster]["phylum"])
            + "\t"
            + str(taxon_paths[cluster]["class"])
            + "\t"
            + str(taxon_paths[cluster]["order"])
            + "\t"
            + str(taxon_paths[cluster]["family"])
            + "\t"
            + str(taxon_paths[cluster]["genus"])
            + "\t"
            + str(taxon_paths[cluster]["species"])
            + "\t"
            + str(top_taxids[cluster])
            + "\n"
        )
    else:
        output_table.write(
            str(cluster)
            + "\t"
            + str(taxon_paths[cluster]["superkingdom"])
            + "\t"
            + str(taxon_paths[cluster]["phylum"])
            + "\t"
            + str(taxon_paths[cluster]["class"])
            + "\t"
            + str(taxon_paths[cluster]["order"])
            + "\t"
            + str(taxon_paths[cluster]["family"])
            + "\t"
            + str(taxon_paths[cluster]["genus"])
            + "\t"
            + str(taxon_paths[cluster]["species"])
            + "\t"
            + str(top_taxids[cluster])
            + "\n"
        )
output_table.close
