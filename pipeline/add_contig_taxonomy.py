#!/usr/bin/env python

# Program that adds contig taxonomy information to a table made by make_contig_table.py
# Uses protein taxonomy information to construct contig taxonomy
# Algorithm:
# In descending order of species, genus, family, order, class, phylum, superkingdom:
#    In descending order of votes:
#        If classification shares common ancestry with majority of other proteins, accept result
#    If no result, move up to next taxonomic level

import sys
from time import gmtime, strftime
import pprint
pp = pprint.PrettyPrinter(indent=4)

rank_priority = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

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
		for current_taxid in contigDictionary[rankName]:
			if isCommonAncestor(current_taxid, taxid, taxidDictionary):
				consistentTotal += 1
			else:
				inconsistentTotal += 1

	if consistentTotal > inconsistentTotal:
		return True
	else:
		return False

def isCommonAncestor(potentialParentTaxid, childTaxid, taxidDictionary):
	current_taxid = childTaxid
	while current_taxid != 1:
		if potentialParentTaxid == current_taxid:
			return True
		current_taxid = taxidDictionary[current_taxid]['parent']
	return False

tax_table_path = sys.argv[1]
taxdump_dir_path = sys.argv[2]
output_file_path = sys.argv[3]

# Process NCBI taxdump files
names_dmp_path = taxdump_dir_path + '/names.dmp'
nodes_dmp_path = taxdump_dir_path + '/nodes.dmp'

taxids = {}
print strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid names'
with open(names_dmp_path) as names_dmp:
	for line in names_dmp:
		line_list = line.rstrip('\n').split('|')
		# Remove trailing and leading spaces
		for i,value in enumerate(line_list):
			line_list[i] = value.strip()

		# Only add scientific name entries
		if line_list[3] == 'scientific name':
			# line_list[1] = line_list[1].replace(' ', '_')
			taxids[line_list[0]] = { 'name': line_list[1] }

print strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid nodes'
with open(nodes_dmp_path) as nodes_dmp:
	for line in nodes_dmp:
		line_list = line.rstrip('\n').split('|')
		# Remove trailing and leading spaces
		for i,value in enumerate(line_list):
			line_list[i] = value.strip()

		taxids[ line_list[0] ][ 'parent' ] = line_list[1]
		taxids[ line_list[0] ][ 'rank' ] = line_list[2]

#pp.pprint(taxids)

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
canonical_ranks = {
	'superkingdom': 1,
	'phylum': 1,
	'class': 1,
	'order': 1,
	'family': 1,
	'genus': 1,
	'species': 1
}
with open(tax_table_path) as tax_table:
	for line in tax_table:
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
			print 'Error, could not find taxid: '
			print line
			sys.exit(2)

		# Now get the taxid of the next canonical rank (if applicable)
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

print strftime("%Y-%m-%d %H:%M:%S") + ' Ranking taxids'
top_taxids = {}
for contig in protein_classifications:
	acceptedTaxid = None
	for rank in rank_priority:
		# Order in descending order of votes
		ordered_taxids = sorted(protein_classifications[contig][rank], key=protein_classifications[contig][rank].__getitem__, reverse=True)
		for taxid in ordered_taxids:
			if isConsistentWithOtherOrfs(taxid, rank, protein_classifications[contig], taxids):
				acceptedTaxid = taxid
				break
	if acceptedTaxid is None:
		acceptedTaxid = 1 # Root

	top_taxids[contig] = acceptedTaxid

pp.pprint(top_taxids)
