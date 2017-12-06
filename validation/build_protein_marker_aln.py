

#input = whole genome fastas
#output = aligned fasta file with all shared protein markers

1. Call ORFs, translate, and identify markers with AMPHORA’s MarkerScanner.pl
2. Combine protein marker sequences with cat
3. Reformat combined protein sequence file with another run of AMPHORA’s MarkerScanner.pl
•This creates .pep files that appear consistent in format with/media/box1/AMPHORA_example/AMPHORA/cobmined/*.pep
4. Align markers with AMPHORA’s MarkerAlignTrim.pl
5. Identify list of shared markers with Python script
6. Concatenate aligned marker sequences with amphora_alignment_filter_and_concat.pl
7. Build tree from alignment with FastTreeMP1. Call ORFs, translate, and identify markers with AMPHORA’s MarkerScanner.pl#wd: /home/ijmiller/AMPHORA2/MIX-51/GCA_000156495.1_ASM15649v1_genomicperl ../Scripts/MarkerScanner.pl -DNA \-Bacteria /home/ijmiller/MIX-51/reference_fasta_files/GCA_000156495.1_ASM15649v1_genomic.fna(Do the same in /home/ijmiller/AMPHORA2/MIX-51//GCA_000012825.1_ASM1282v1_genomic2. Combine protein marker sequences with cat#wd: /home/ijmiller/AMPHORA2/MIX-51/combinedcat ../GCA_000156495.1_ASM15649v1_genomic/*pep \../GCA_000012825.1_ASM1282v1_genomic/*pep>combined.pep3. Reformat combined protein sequence file with another run of AMPHORA’s MarkerScan-ner.plperl ../../Scripts/MarkerScanner.pl -Bacteria combined.pepThis effectivelly separates out all of the markers into individual multifasta .pep files that appear consistent informat with media/box1/AMPHORA_example/AMPHORA/cobmined/*.pep1
#1. Run


#4. Concatenate .aln files

#Testing in: /home/ijmiller/AMPHORA2/MIX-51/combined
#dnaG.aln has 3 sequences and rplT.aln has 1

#4a. Identify lift of shared sequences

!rm combined.aln

import glob
import re

num_input_genome = 2
shared_protein_list = []
genome_dict = {}
num_marker_proteins = 0
for aln in glob.glob("*.aln"):
    num_marker_proteins += 1
    protein_name = aln.rstrip(".aln")
    with open(aln) as infile:
        num_aln_genomes = 0
        for line in infile:
            if ">" in line:
                seq_name = line.rstrip().lstrip(">")
                genome = seq_name.split("SENSE) ")[-1]
                m = re.search(r"[A-Z]{1}[a-z]+\s[a-z]+", seq_name)
                if m:
                    genome = m.group()
                else:
                    print("No match found:\n{}".format(seq_name))

                if genome not in genome_dict:
                    genome_dict[genome] = {}
                    genome_dict[genome][protein_name] = ""
                    num_aln_genomes += 1
                elif protein_name not in genome_dict[genome]:
                    genome_dict[genome][protein_name] = ""
                    num_aln_genomes += 1
                else:
                    print("Duplicate protein for {}\nOmitting duplicate...".format(protein_name))

            else:
                genome_dict[genome][protein_name] += line.rstrip()
        if num_aln_genomes == num_input_genome:
            #This will not properly handle the case of multiple protein
            #seqs from single genome
            shared_protein_list.append(protein_name)
        else:
            print("{} not presented in all {} genomes, omitting from tree...".format(protein_name,num_input_genome))

print("\nUsing {} of {} protein markers to build concatenated alignment...".format(num_marker_proteins,len(shared_protein_list)))
combined_aln_dict = {}
for genome in genome_dict.keys():
    combined_aln_dict[genome] = ""
    for count,shared_protein in enumerate(shared_protein_list):
        print genome,shared_protein
        seq = genome_dict[genome][shared_protein]
        corrected_seq = ""
        #Replace masked letters
        for letter in seq:
            if letter.islower():
                corrected_seq += "-"
            else:
                corrected_seq += letter
        combined_aln_dict[genome] += "-" * 10 + corrected_seq

with open("combined.aln", "w") as outfile:
    for genome in combined_aln_dict:
        outfile.write(">{}\n".format(genome))
        outfile.write("{}\n".format(combined_aln_dict[genome]))
