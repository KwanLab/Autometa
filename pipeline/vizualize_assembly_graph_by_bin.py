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

# This program uses a spades assembly graph and a binning table to produce
# files to visualize the graph in cytoscape, with contigs colored by bin

import argparse
import os
import pandas as pd
import itertools
import gzip
import pdb


def getGraph(graph_file, paths_file):
    # Load graph file into memory
    S_lines = list()
    L_lines = list()
    P_lines = list()

    if os.path.splitext(graph_file)[-1] == ".gz":
        input_graph = gzip.open(graph_file, "rb")
    else:
        input_graph = open(graph_file)

    for line in input_graph:
        line_list = line.rstrip().split("\t")
        if line_list[0] == "S":
            S_lines.append(line)
        if line_list[0] == "L":
            L_lines.append(line)
        if line_list[0] == "P":
            P_lines.append(line)

    input_graph.close()

    # First we get the lengths of each segment, which will be useful later
    segment_lengths = dict()  # Keyed by segment(int)
    for line in S_lines:
        line_list = line.rstrip().split("\t")
        segment_name = int(line_list[1])
        segment_length = len(line_list[2])
        segment_lengths[segment_name] = segment_length

    # We need to identify the end pieces of each scaffold
    scaffolds = (
        list()
    )  # Holds all scaffold names encountered (so that we can make the internal connections between s and e)
    end_segments = (
        dict()
    )  # Keyed by segment + s|e, holds list of scaffold names in the form scaffold_name + s|e
    # Note: We count as an end anything within 100 bp of the end of the scaffold, so a scaffold can have multiple ends

    scaffold_paths = (
        dict()
    )  # Keyed by scaffold name, holds lists of tuples in the form (segment(int), orientation (+|-))
    # Here we diverge depending on whether there is a seperate paths file or not
    if paths_file:
        current_scaffold = None
        record = False
        current_seg_list = list()

        if os.path.splitext(paths_file)[-1] == ".gz":
            paths = gzip.open(paths_file, "rb")
        else:
            paths = open(paths_file)

        for line in paths:
            if len(line) > 5 and line[0:5] == "NODE_" and line.rstrip()[-1] != "'":
                # First record the last list
                scaffold_paths[current_scaffold] = current_seg_list

                # This is the start of a scaffold
                record = True
                current_scaffold = line.rstrip()
                scaffolds.append(current_scaffold)
                current_seg_list = list()
            elif len(line) > 5 and line[0:5] == "NODE_" and line.rstrip()[-1] == "'":
                record = False
            else:
                if record == True:
                    clean_line = line.rstrip(";\n")
                    initial_path_list = clean_line.split(",")
                    for segment_string in initial_path_list:
                        segment = int(segment_string[:-1])
                        orientation = segment_string[-1]
                        segment_tuple = (segment, orientation)
                        current_seg_list.append(segment_tuple)
        # Record last path list
        scaffold_paths[current_scaffold] = current_seg_list

        paths.close()
    else:
        # In this case we get the paths directly from the gfa file
        for line in P_lines:
            line_list = line.rstrip().split("\t")
            scaffold_name = "_".join(line_list[1].split("_")[:-1])
            scaffolds.add(scaffold_name)
            initial_path_list = line_list[2].split(",")
            scaffold_path_list = list()
            for segment_string in initial_path_list:
                segment = int(segment_string[:-1])
                orientation = segment_string[-1]
                segment_tuple = (segment, orientation)
                scaffold_path_list.append(segment_tuple)
            if scaffold_name in scaffold_paths:
                scaffold_paths[scaffold_name].extend(scaffold_path_list)
            else:
                scaffold_paths[scaffold_name] = scaffold_path_list

    # Now we work out which are the end segments
    for scaffold_name in scaffold_paths:
        scaffold_path_list = scaffold_paths[scaffold_name]
        # We can assume here that the list has at least one member, but we cannot assume it has at least two

        seg_ends = list()
        scaffold_ends = list()

        # If the end segment on either side is < 100 bp then we add some more end tuples
        if len(scaffold_path_list) > 1:
            length_traversed = 0
            for i in range(0, len(scaffold_path_list)):
                if length_traversed < 100:
                    if scaffold_path_list[1][1] == "+":
                        seg_ends.append(str(scaffold_path_list[i][0]) + "s")
                    else:
                        seg_ends.append(str(scaffold_path_list[i][0]) + "e")
                    scaffold_ends.append(scaffold_name + "s")
                else:
                    break
                try:
                    length_traversed += segment_lengths[scaffold_path_list[i][0]]
                except:
                    pdb.set_trace()

            length_traversed = 0
            for i in reversed(list(range(0, len(scaffold_path_list)))):
                if length_traversed < 100:
                    if segment_lengths[scaffold_path_list[i][1]] == "+":
                        seg_ends.append(str(scaffold_path_list[i][0]) + "e")
                    else:
                        seg_ends.append(str(scaffold_path_list[i][0]) + "s")
                    scaffold_ends.append(scaffold_name + "e")
                else:
                    break

                length_traversed += segment_lengths[scaffold_path_list[i][0]]

        for i in range(len(seg_ends)):
            if seg_ends[i] in end_segments:
                end_segments[seg_ends[i]].append(scaffold_ends[i])
            else:
                end_segments[seg_ends[i]] = [scaffold_ends[i]]

    # Now we start making a represenation of the graph in memory

    # First, intra connections (between s and e of same scaffold)
    graph = dict()  # Dictionary of lists
    for scaffold in scaffolds:
        graph[scaffold + "s"] = [scaffold + "e"]
        graph[scaffold + "e"] = [scaffold + "s"]

    # Now we make connections on the basis of two contig ends sharing the same segment
    end_segment_names = set()
    for segment_name in end_segments:
        end_segment_names.add(segment_name[:-1])  # i.e. cut off end 's' or 'e'

    shared_segments = set()
    for segment_name in end_segment_names:
        if segment_name + "s" in end_segments and segment_name + "e" in end_segments:
            shared_segments.add(segment_name)

    for segment in shared_segments:
        node_list1 = end_segments[segment + "s"]
        node_list2 = end_segments[segment + "e"]

        # See logic's answer to stack overflow query at https://stackoverflow.com/questions/12935194/combinations-between-two-lists
        combination_tuple_list = [(x, y) for x in node_list1 for y in node_list2]

        for combination_tuple in combination_tuple_list:
            # First check whether link exists
            if combination_tuple[1] not in graph[combination_tuple[0]]:
                graph[combination_tuple[0]].append(combination_tuple[1])

            if combination_tuple[0] not in graph[combination_tuple[1]]:
                graph[combination_tuple[1]].append(combination_tuple[0])

    # Now we add the inter connections from the L lines of the gfa file
    for line in L_lines:
        line_list = line.rstrip().split("\t")
        seg1 = line_list[1]
        seg1_orientation = line_list[2]
        seg2 = line_list[3]
        seg2_orientation = line_list[4]

        if seg1_orientation == "+":
            segment1_end = seg1 + "e"
        else:
            segment1_end = seg1 + "s"

        if seg2_orientation == "+":
            segment2_end = seg2 + "s"
        else:
            segment2_end = seg2 + "e"

        # First check that both are end segments
        if segment1_end in end_segments and segment2_end in end_segments:
            # We make a list of nodes that are associated with both segment1_end and segment2_end, and output lines for each combination
            node_list1 = end_segments[segment1_end]
            node_list2 = end_segments[segment2_end]

            # See above
            combination_tuple_list = [(x, y) for x in node_list1 for y in node_list2]

            for combination_tuple in combination_tuple_list:
                # First check whether link exists
                if combination_tuple[1] not in graph[combination_tuple[0]]:
                    graph[combination_tuple[0]].append(combination_tuple[1])

                if combination_tuple[0] not in graph[combination_tuple[1]]:
                    graph[combination_tuple[1]].append(combination_tuple[0])

    return graph


def bfs(graph, start_set):
    # Expects graph to be a dictionary of dictionaries, and start to be a set of starting nodes
    # The start set has bare contig names, whereas the graph has contig + 's'|'e'

    # Keep track of all visited nodes
    explored = []
    # Keep track of nodes to be checked
    queue = []
    for contig in start:
        queue.append(contig + "s")
        queue.append(contig + "e")

    # Keep looping until there are no nodes still to be checked
    while queue:
        # pop shallowest node (first node) from queue
        node = queue.pop(0)
        if node not in explored:
            # add node to list of checked nodes
            explored.append(node)
            if node in graph:
                for neighbor in graph[node]:
                    if neighbor not in explored:
                        queue.append(neighbor)

    return set(explored)


parser = argparse.ArgumentParser(
    description="Script to parse a SPAdes assembly graph in order to visualize a binned metagenome in cytoscape, allowing bin refinement. Uses breadth-first search to make sub-files with less complexity."
)
parser.add_argument(
    "-b",
    "--bin_table",
    metavar="<bin.tab>",
    help="path to the output from either run_autometa.py or ML_recruitment.py",
    required=True,
)
parser.add_argument(
    "-c",
    "--column",
    metavar="<bin column name>",
    help="the name of the column to use for binning purposes",
    default="cluster",
)
parser.add_argument(
    "-o",
    "--output_dir",
    metavar="<dir>",
    help="path to the directory where output files will go",
    default=".",
)
parser.add_argument(
    "-g",
    "--graph_file",
    metavar="<assembly_graph_with_scaffolds.gfa|assembly_graph.gfa>",
    help="Path to the assembly graph from SPAdes. Use either assembly_graph_with_scaffolds.gfa or assembly_graph.gfa (if you used an old version of SPAdes that doesn't make an assembly_graph_with_scaffolds.gfa file). If the latter, you MUST provide a scaffolds.paths file",
    required=True,
)
parser.add_argument(
    "-p",
    "--paths_file",
    metavar="<scaffolds.paths>",
    help="Path to the scaffolds.paths file made by SPAdes",
)

args = vars(parser.parse_args())

bin_table_path = os.path.abspath(args["bin_table"])
cluster_column_heading = args["column"]
output_dir = os.path.abspath(args["output_dir"])
graph_file_path = os.path.abspath(args["graph_file"])
paths_file_path = os.path.abspath(args["paths_file"])

# Check paths exist
if not os.path.isfile(bin_table_path):
    print(
        ("Error! Could not find a bin table at the following path: " + bin_table_path)
    )
    exit(1)

if not os.path.isfile(graph_file_path):
    print(("Error! Cannot find a graph file at the following path: " + graph_file_path))
    exit(1)

if paths_file_path and not os.path.isfile(paths_file_path):
    print(("Error! Cannot find a paths file at the following path: " + paths_file_path))
    exit(1)

graph_filename = graph_file_path.split("/")[-1]
if not (
    graph_filename == "assembly_graph.gfa" or graph_filename == "assembly_graph.gfa.gz"
) and not (
    graph_filename == "assembly_graph_with_scaffolds.gfa"
    or graph_filename == "assembly_graph_with_scaffolds.gfa.gz"
):
    print(
        "Error! You must provide either the file assembly_graph.gfa or assembly_graph_with_scaffolds.gfa as the graph file"
    )
    exit(1)

# Check that a paths file is supplied if the graph file is assembly_graph.gfa
if (
    graph_filename == "assembly_graph.gfa" or graph_filename == "assembly_graph.gfa.gz"
) and paths_file_path is None:
    print(
        "Error! If you provide the file assembly_graph.gfa, you must also provide scaffolds.paths"
    )
    exit(1)

# Make output directory if it doesn't exist
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

master_table = pd.read_table(bin_table_path)

# Format check for the table
columns_to_check = [cluster_column_heading, "contig"]

for column in columns_to_check:
    if column not in master_table.columns:
        print(
            (
                "Error! Could not find a column called "
                + column
                + " in table "
                + bin_table_path
            )
        )
        exit(1)

# Here is what we are going to do:
# Make a datastructure of the contigs in each bin
# Output a table that can be used to color nodes in cytoscape
# Make a representation of the assembly graph in memory
# Output the full graph file
# Output a subgraph for all bins, identified through BFS from each bin
# We check the BFS results to see if there is overlap in bins. If there is, we output one file for n bins that overlap

bin_contigs = dict()  # keyed by bin name, holds lists of contigs
for i, row in master_table.iterrows():
    contig = row["contig"]
    bin_name = row[cluster_column_heading]
    if bin_name == "unclustered":
        continue
    if bin_name in bin_contigs:
        bin_contigs[bin_name].append(contig)
    else:
        bin_contigs[bin_name] = [contig]

color_table_path = os.path.join(output_dir, "color_table")
color_table = open(color_table_path, "w")
color_table.write("contig\tbin\n")

for bin_name in bin_contigs:
    for contig in bin_contigs[bin_name]:
        color_table.write(contig + "s\t" + bin_name + "\n")
        color_table.write(contig + "e\t" + bin_name + "\n")
color_table.close()

# Load the graph
assembly_graph = getGraph(graph_file_path, paths_file_path)

# We now make subgraphs for each bin (by BFS)
bin_bfs_sets = dict()  # Keyed by bin, holds BFS node sets
for bin_name in bin_contigs:
    bin_bfs_sets[bin_name] = bfs(assembly_graph, bin_contigs[bin_name])

# Now we work out which BFS sets overlap
merge_done = True
new_bfs_sets = dict()  # Will hold new bfs_sets
while list(bin_bfs_sets.keys()):
    bin_list = list(bin_bfs_sets.keys())
    bin_to_check = bin_list.pop(0)
    bins_to_merge = []
    for bin_name in bin_list:
        if bool(bin_bfs_sets[bin_to_check] & bin_bfs_sets[bin_name]):
            bins_to_merge.append(bin_name)
    # Merge bins if needed
    if len(bins_to_merge) > 0:
        bins_to_merge = [bin_to_check] + bins_to_merge
        new_bin_name = ",".join(bins_to_merge)
        new_bin_set = set()
        for bin_name in bins_to_merge:
            new_bin_set = new_bin_set.union(bin_bfs_sets[bin_name])
        # Delete old bfs sets
        for bin_name in bins_to_merge:
            del bin_bfs_sets[bin_name]
        # Make new bin
        bin_bfs_sets[new_bin_name] = new_bin_set
    else:
        # Move bin to the new data structure
        new_bfs_sets[bin_to_check] = bin_bfs_sets[bin_to_check]
        del bin_bfs_sets[bin_to_check]

# Now we output graphs for each merged bfs set
filehandles = dict()  # Keyed by new bin name
nodes_explored = (
    set()
)  # Protects from writing a link if both nodes have already been explored
# (i.e. cuts down on file redundancy)

for bin_name in new_bfs_sets:
    output_path = os.path.join(output_dir, bin_name + "_graph.txt")
    filehandles[bin_name] = open(output_path, "w")
    filehandles[bin_name].write("node1\ttype\tnode2\tedge\n")

for source_node in assembly_graph:
    destination_node_list = assembly_graph[source_node]
    for bin_name in new_bfs_sets:
        if source_node in new_bfs_sets[bin_name]:
            for destination_node in destination_node_list:
                if (
                    source_node not in nodes_explored
                    and destination_node not in nodes_explored
                ):
                    # Write to relevant file
                    if source_node[:-1] == destination_node[:-1]:
                        edge_type = "intra"
                    else:
                        edge_type = "inter"
                    filehandles[bin_name].write(
                        source_node
                        + "\t"
                        + edge_type
                        + "\t"
                        + destination_node
                        + "\t1\n"
                    )
                    nodes_explored.add(source_node)
                    nodes_explored.add(destination_node)

# Close files
for bin_name in filehandles:
    filehandles[bin_name].close()
