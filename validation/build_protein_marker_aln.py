#!/usr/bin/env python

import glob
import os
import re
import subprocess
import pandas as pd
from collections import Counter
import argparse

parser = argparse.ArgumentParser(description='Script to build an aligned fasta\
    file using AMPHORA2 from a directory of fasta files.')
parser.add_argument('-x','--extension', help='File extension to use for fasta files', default=".fasta")
args = vars(parser.parse_args())

#input = whole genome fastas
#output = aligned fasta file with all shared protein markers
"""
0. Standardize fasta contig names
1. Call ORFs, translate, and identify markers with AMPHORA's MarkerScanner.pl
2. Combine protein marker sequences with cat
3. Reformat combined protein sequence file with another run of AMPHORA's MarkerScanner.pl
    This creates .pep files that appear consistent in format with/media/box1/AMPHORA_example/AMPHORA/cobmined/*.pep
4. Align markers with AMPHORA's MarkerAlignTrim.pl
5. Identify list of shared markers and extract individual aligned sequences
6. Correct and concatenate aligned marker sequences
7. Build tree from alignment with FastTreeMP.
"""

for fna in glob.glob("*" + args['extension']):
    basename = ".".join(fna.split(".")[:-1])
    output_handle = basename + "_renamed.fasta"
    input_handle = open(fna)
    with open(output_handle,"w") as outfile:
        for count,line in enumerate(input_handle):
            if ">" in line:
                seq_name = line.rstrip().lstrip(">")
                outfile.write(">" + basename + "-" + str(count) + "\n")
            else:
                outfile.write(line)
    input_handle.close()

#1. For every fasta in the working directory
num_input_genome = 0
for fasta in glob.glob("*_renamed.fasta"):
    #Make a new directory with the fastas base name, if it doesn't exist
    basename = ".".join(fasta.split(".")[:-1])
    if not os.path.isdir(basename):
        subprocess.call(["mkdir",basename])
    #Run MarkerScanner.pl in the new folder
    os.chdir(basename)
    subprocess.call(["MarkerScanner.pl", "-DNA","-Bacteria","../"+fasta])
    os.chdir("../")
    num_input_genome += 1

#2. Combine protein marker sequences with cat
#subprocess.call("find . -mindepth 2 -name '*.pep'  | xargs -I % cat % >> combined.pep", shell=True)
subprocess.call("find . -mindepth 2 -maxdepth 2 -name '*.pep'  | xargs -I % cat % >> combined.pep", shell=True)

#3. Reformat combined protein sequence file with another run of AMPHORA's MarkerScanner.pl
subprocess.call(["mkdir","combined"])
os.chdir("combined") #/home/ijmiller/AMPHORA2/MIX-51/fasta_test/combined
subprocess.call(["MarkerScanner.pl","-Bacteria","../combined.pep"])

#4. Align markers with AMPHORA's MarkerAlignTrim.pl
subprocess.call(["MarkerAlignTrim.pl","-OutputFormat","fasta"])
aln_list = glob.glob("*.aln")
if os.path.exists("combined.aln"):
    aln_list.remove("combined.aln")

#5.  Identify list of shared markers and extract individual aligned sequences
#num_input_genome = 52 #--> Defined above
bad_contig_list = []
shared_protein_list = []
marker_list = []
genome_dict = {}
num_marker_proteins = 0
skip_sequence = False
for aln in aln_list:
    num_marker_proteins += 1
    protein_name = aln.rstrip(".aln")
    marker_list.append(protein_name)
    with open(aln) as infile:
        num_aln_genomes = 0
        aln_dict = {}
        for line in infile:

            if ">" in line:
                skip_sequence = False
                seq_name = line.rstrip().lstrip(">")
                genome = seq_name.split("SENSE) ")[-1]
                #m = re.search(r">[A-Z]+.*?-[0-9]",line.rstrip())
                #m = re.search(r">cluster_DBSCAN_round.*?-",line.rstrip())
                m = re.search(r">.*?-",line.rstrip())
                if m:
                    genome = m.group()[1:-1]
                    #genome = genome_filename_dict[genome_basename]
                else:
                    print("No match found:\n{}\nExitting...".format(seq_name))
                    print("Make sure to delete files before rerunning AMPHORA scripts...")
                    bad_contig_list.append(seq_name)
                    break

                if genome not in genome_dict:
                    genome_dict[genome] = {}
                    genome_dict[genome][protein_name] = ""
                    num_aln_genomes += 1
                    aln_dict[protein_name]=[genome]
                elif protein_name not in genome_dict[genome]:
                    genome_dict[genome][protein_name] = ""
                    num_aln_genomes += 1
                    #aln_dict[protein_name].append(genome)
                else:
                    skip_sequence = True
                    print("Ommitting duplicate {} protein for {}...".format(protein_name,genome))

            elif skip_sequence: #Skip duplicate markers
                continue
            else:
                genome_dict[genome][protein_name] += line.rstrip()

        if num_aln_genomes == num_input_genome:
        #if len(aln_dict.keys()) == num_input_genome:
            #This will not properly handle the case of multiple protein
            #seqs from single genome
            shared_protein_list.append(protein_name)

        else:
            #print("{} genomes have {}".format(num_aln_genomes,protein_name))
            print("{} not presented in all {} genomes, omitting from tree...".format(protein_name,num_input_genome))
print("\nUsing {} of {} protein markers to build concatenated alignment...".format(len(shared_protein_list),num_marker_proteins))

if len(shared_protein_list) == 0:
    for genome in genome_dict.keys():
        print genome,len(genome_dict[genome].keys())
        missing_markers = [marker for marker in marker_list if marker not in genome_dict[genome].keys()]
        if len(missing_markers) > 0:
            print ("\tMissing: {}".format(",".join(missing_markers)))

#NOTE: Should check length - Jason's note in 'amphora_merge_alignments.pl'
# Sometimes APMPHORA2 misses a character on the end of the last sequence
# Check the alignment lengths
print("\nChecking to make sure all shared markers have same alignment length:")
aln_length_dict = {}
for shared_protein in shared_protein_list:
    for genome in genome_dict.keys():
        protein_aln_len = len(genome_dict[genome][shared_protein])
        if not shared_protein in aln_length_dict:
            aln_length_dict[shared_protein] = [protein_aln_len]
        else:
            aln_length_dict[shared_protein].append(protein_aln_len)
    print shared_protein,Counter(aln_length_dict[shared_protein])

#6. Correct and concatenate aligned marker sequences
combined_aln_dict = {}
for genome in genome_dict.keys():
    combined_aln_dict[genome] = ""
    for count,shared_protein in enumerate(shared_protein_list):
        #print genome,shared_protein
        seq = genome_dict[genome][shared_protein]
        corrected_seq = ""
        #Replace masked letters
        for letter in seq:
            if letter.islower():
                corrected_seq += "-"
            #Replace '*' with '-'
            elif letter == "*":
                corrected_seq += "-"
            else:
                corrected_seq += letter

        combined_aln_dict[genome] += "-" * 10 + corrected_seq

with open("combined.aln", "w") as outfile:
    for genome in combined_aln_dict:
        outfile.write(">{}\n".format(genome))
        outfile.write("{}\n".format(combined_aln_dict[genome]))

#7. Build tree from alignment with FastTreeMP.
