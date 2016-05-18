#!/usr/bin/env python

# Program that adds contig taxonomy information to a table made by make_contig_table.py
# Uses protein taxonomy information to construct contig taxonomy
# Algorithm:
# In descending order of species, genus, family, order, class, phylum, superkingdom:
#    In descending order of votes:
#        If classification shares common ancestry with majority of other proteins, accept result
#    If no result, move up to next taxonomic level

import sys
from time import *
from tqdm import *
import subprocess
import pprint
import pdb
pp = pprint.PrettyPrinter(indent=4)

rank_priority = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'root']
canonical_ranks = {
	'superkingdom': 1,
	'phylum': 1,
	'class': 1,
	'order': 1,
	'family': 1,
	'genus': 1,
	'species': 1
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
		current_taxid = taxidDictionary[current_taxid]['parent']
	return False

def lowest_majority(contigDictionary, taxidDictionary):
	taxid_totals = {} # Dictionary of dictionary, keyed by rank then by taxid, holds totals accounting for whole taxid paths
	for rank in rank_priority:
		if rank in contigDictionary:
			for taxid in contigDictionary[rank]:
				# We need to add to taxid_totals for each taxid in the tax_path
				current_taxid = taxid
				current_rank = rank
				while int(current_taxid) != 1:
					if current_rank in canonical_ranks:
						if current_rank in taxid_totals:
							if current_taxid in taxid_totals[current_rank]:
								taxid_totals[current_rank][current_taxid] += 1
							else:
								taxid_totals[current_rank][current_taxid] = 1
						else:
							taxid_totals[current_rank] = { current_taxid: 1 }
					current_taxid = taxidDictionary[current_taxid]['parent']
					current_rank = taxidDictionary[current_taxid]['rank']

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
			majority_threshold = float(total_votes)/2
			if taxid_leader_votes > majority_threshold:
				return taxid_leader_votes

	# Just in case
	return 1

contig_table_path = sys.argv[1]
tax_table_path = sys.argv[2]
taxdump_dir_path = sys.argv[3]
output_file_path = sys.argv[4]

# Process NCBI taxdump files
names_dmp_path = taxdump_dir_path + '/names.dmp'
nodes_dmp_path = taxdump_dir_path + '/nodes.dmp'

taxids = {}
print strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid names'
wc_output = subprocess.check_output(['wc', '-l', names_dmp_path])
wc_list = wc_output.split()
number_of_lines = int(wc_list[0])

with open(names_dmp_path) as names_dmp:
	for line in tqdm(names_dmp, total=number_of_lines):
		line_list = line.rstrip('\n').split('|')
		# Remove trailing and leading spaces
		for i,value in enumerate(line_list):
			line_list[i] = value.strip()

		# Only add scientific name entries
		if line_list[3] == 'scientific name':
			# line_list[1] = line_list[1].replace(' ', '_')
			taxids[line_list[0]] = { 'name': line_list[1] }

print strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid nodes'
wc_output = subprocess.check_output(['wc', '-l', nodes_dmp_path])
wc_list = wc_output.split()
number_of_lines = int(wc_list[0])

with open(nodes_dmp_path) as nodes_dmp:
	for line in tqdm(nodes_dmp, total=number_of_lines):
		line_list = line.rstrip('\n').split('|')
		# Remove trailing and leading spaces
		for i,value in enumerate(line_list):
			line_list[i] = value.strip()

		taxids[ line_list[0] ][ 'parent' ] = line_list[1]
		taxids[ line_list[0] ][ 'rank' ] = line_list[2]

name_lookup = {} # Dictionary of dictionaries, keyed by rank then name
for taxid in taxids:
	rank = taxids[taxid]['rank']
	name = taxids[taxid]['name']
	if rank in name_lookup:
		name_lookup[rank][name] = taxid
	else:
		name_lookup[rank] = { name: taxid }

print strftime("%Y-%m-%d %H:%M:%S") + ' Parsing taxonomy table'
protein_classifications = {} # protein_classifications[contig][rank][taxid] (running total of each thing)
number_of_proteins = {}

# Work out number of lines in file
wc_output = subprocess.check_output(['wc', '-l', tax_table_path])
wc_list = wc_output.split()
number_of_lines = int(wc_list[0])

with open(tax_table_path) as tax_table:
	for line in tqdm(tax_table, total=number_of_lines):

		line_list = line.rstrip('\n').split('\t')
		seqNameList = line_list[0].split('_')
		seqNameList.pop()
		contigName = ('_').join(seqNameList)

		# Get taxid
		taxName = line_list[1]
		taxRank = line_list[2]
		taxid = None
		if taxRank in name_lookup.keys() and taxName in name_lookup[taxRank].keys():
			taxid = name_lookup[taxRank][taxName]
		else:
			taxid = 1 # Treat unknowns as root

		# Now get the taxid of the next canonical rank (if applicable)
		if taxRank == 'no rank':
			taxid = 1
			taxRank = 'root'

		if taxid != 1:
			while taxRank not in rank_priority:
				taxid = taxids[taxid]['parent']
				taxRank = taxids[taxid]['rank']

		# Keep running total of taxids for each contig
		if contigName not in protein_classifications:
			protein_classifications[contigName] = {}
		if taxRank not in protein_classifications[contigName]:
			protein_classifications[contigName][taxRank] = {}

		if taxid not in protein_classifications[contigName][taxRank]:
			protein_classifications[contigName][taxRank][taxid] = 1
		else:
			protein_classifications[contigName][taxRank][taxid] += 1

		# Count number of proteins per contig
		if contigName in number_of_proteins:
			number_of_proteins[contigName] += 1
		else:
			number_of_proteins[contigName] = 1

pdb.set_trace()

print strftime("%Y-%m-%d %H:%M:%S") + ' Ranking taxids'
top_taxids = {}
total_contigs = len(protein_classifications)

for contig in tqdm(protein_classifications, total=total_contigs):
	acceptedTaxid = None
	for rank in rank_priority:
		if acceptedTaxid is not None:
			break
		# Order in descending order of votes
		if rank in protein_classifications[contig]:
			ordered_taxids = sorted(protein_classifications[contig][rank], key=protein_classifications[contig][rank].__getitem__, reverse=True)
			#sys.exit()
			for taxid in ordered_taxids:
				if isConsistentWithOtherOrfs(taxid, rank, protein_classifications[contig], taxids):
					acceptedTaxid = taxid
					break
	
	# If acceptedTaxid is still None at this point, there was some kind of draw, so we need to find the lowest taxonomic level where there is a
	# majority
	if acceptedTaxid is None:
		acceptedTaxid = lowest_majority(protein_classifications[contig], taxids)

	top_taxids[contig] = acceptedTaxid

pdb.set_trace()

print strftime("%Y-%m-%d %H:%M:%S") + ' Resolving taxon paths'
taxon_paths = {} # Dictionary of dictionaries, keyed by contig then rank, contains the taxon names
for contig in tqdm(top_taxids, total=total_contigs):
	taxon_paths[contig] = {}
	current_taxid = top_taxids[contig]

	while current_taxid != 1:
		current_rank = taxids[current_taxid]['rank']
		if current_rank in canonical_ranks:
			taxon_paths[contig][current_rank] = taxids[current_taxid]['name']
		current_taxid = taxids[current_taxid]['parent']

	for rank in rank_priority:
		if rank not in taxon_paths[contig]:
			taxon_paths[contig][rank] = 'unclassified'

pdb.set_trace()

print strftime("%Y-%m-%d %H:%M:%S") + ' Writing table'
output_table = open(output_file_path, 'w')
with open(contig_table_path) as contig_table:
	for i,line in enumerate(tqdm(contig_table, total=(total_contigs+1))):
		if i == 0:
			original_line = line.rstrip('\n')
			new_header = original_line + '\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\ttaxid\n'
			output_table.write(new_header)
		else:
			original_line = line.rstrip('\n')
			line_list = original_line.split('\t')
			contig_name = line_list[0]
			if contig_name not in top_taxids:
				print 'Error, could not find ' + contig_name + ' in ' + output_file_path
				sys.exit(2)
			new_line = origina_line + '\t' + taxon_paths[contig_name]['kingdom'] + '\t' + taxon_paths[contig_name]['phylum'] + '\t' + taxon_paths[contig_name]['class'] + '\t' + taxon_paths[contig_name]['order'] + '\t' + taxon_paths[contig_name]['family'] + '\t' + taxon_paths[contig_name]['genus'] + taxon_paths[contig_name]['species'] + '\t' + top_taxids[contig_name] + '\n'
			output_table.write(new_line)
output_table.close
