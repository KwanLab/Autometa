#!/usr/bin/env python3

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

# Program that adds contig taxonomy information to a table made by make_contig_table.py
# Uses protein taxonomy information to construct contig taxonomy
# Algorithm:
# In descending order of species, genus, family, order, class, phylum, superkingdom:
#    In descending order of votes:
#        If classification shares common ancestry with majority of other proteins, accept result
#    If no result, move up to next taxonomic level
# Uses output from program blast2lca (https://github.com/emepyc/Blast2lca)
# Updated version uses output from lca.py in Autometa pipeline


import os
import sys
import subprocess
import time

from tqdm import tqdm

CANONICAL_RANKS = [
    'species',
    'genus',
    'family',
    'order',
    'class',
    'phylum',
    'superkingdom',
    'root'
]

def isConsistentWithOtherOrfs(taxid, rank, ctg_lcas, nodes_dict):
    """
    Function that determines for a given taxid, whether the majority of proteins
    in a contig, with rank equal to or above the given rank, are common
    ancestors of the taxid.  If the majority are, this function returns True,
    otherwise it returns False
    """
    # ctg_lcas = {ctg:{canonical_rank:{taxid:num_hits,...},...},ctg2:{...},...}
    # First we make a modified CANONICAL_RANKS list that only includes the current rank and above
    rank_index = CANONICAL_RANKS.index(rank)
    ranks_to_consider = CANONICAL_RANKS[rank_index:]

    # Now we total up the consistent and inconsistent ORFs
    consistent = 0
    inconsistent = 0

    for rankName in ranks_to_consider:
        if rankName not in ctg_lcas:
            continue
        for ctg_lca in ctg_lcas[rankName]:
            if isCommonAncestor(ctg_lca, taxid, nodes_dict):
                consistent += ctg_lcas[rankName][ctg_lca]
            else:
                inconsistent += ctg_lcas[rankName][ctg_lca]

    if consistent > inconsistent:
        return True
    else:
        return False

def isCommonAncestor(parent_taxid, child_taxid, nodes_dict):
    ancestor_taxid = child_taxid
    while ancestor_taxid != 1:
        if parent_taxid == ancestor_taxid:
            return True
        ancestor_taxid = nodes_dict[ancestor_taxid]['parent']
    return False

def lowest_majority(ctg_lcas, nodes_dict):
    # ctg_lcas = {canonical_rank:{taxid:num_hits, taxid2:#,...},rank2:{...},...}
    taxid_totals = {}

    for rank in CANONICAL_RANKS:
        if rank not in ctg_lcas:
            continue

        rank_index = CANONICAL_RANKS.index(rank)
        ranks_to_consider = CANONICAL_RANKS[rank_index:]

        for taxid in ctg_lcas[rank]:
            # Make a dictionary to total the number of canonical ranks hit
            # while traversing the path so that we can add 'unclassified' to
            # any that don't exist. Later we need to make sure that
            # 'unclassified' doesn't ever win
            ranks_in_path = {rank_to_consider:0 for rank_to_consider in ranks_to_consider}

            # We need to add to taxid_totals for each taxid in the tax_path
            current_taxid = taxid
            current_rank = rank
            while current_taxid != 1:
                if current_rank not in set(CANONICAL_RANKS):
                    current_taxid = nodes_dict[current_taxid]['parent']
                    current_rank = nodes_dict[current_taxid]['rank']
                    continue

                ranks_in_path[current_rank] += 1

                if current_rank not in taxid_totals:
                    taxid_totals.update({current_rank:{current_taxid:1}})
                    current_taxid = nodes_dict[current_taxid]['parent']
                    current_rank = nodes_dict[current_taxid]['rank']
                    continue

                if current_taxid in taxid_totals[current_rank]:
                    taxid_totals[current_rank][current_taxid] += 1
                else:
                    taxid_totals[current_rank][current_taxid] = 1

                current_taxid = nodes_dict[current_taxid]['parent']
                current_rank = nodes_dict[current_taxid]['rank']

            # Now go through ranks_in_path. Where total = 0, add 'unclassified'
            for rank_to_consider in ranks_to_consider:
                if ranks_in_path[rank_to_consider] == 0:
                    if rank_to_consider not in taxid_totals:
                        taxid_totals[rank_to_consider] = {'unclassified':1}
                        continue
                    if 'unclassified' in taxid_totals[rank_to_consider]:
                        taxid_totals[rank_to_consider]['unclassified'] += 1
                    else:
                        taxid_totals[rank_to_consider]['unclassified'] = 1


    # If there are any gaps in the taxonomy paths for any of the proteins in the contig,
    # we need to add 'unclassified' to the relevant canonical taxonomic rank.
    # However, we must never allow 'unclassified' to win! (That just won't really tell us anything)

    # Now we need to determine which is the first level to have a majority
    for rank in CANONICAL_RANKS:
        total_votes = 0
        taxid_leader = None
        taxid_leader_votes = 0
        if not rank in taxid_totals:
            continue
        for taxid in taxid_totals[rank]:
            taxid_votes = taxid_totals[rank][taxid]
            total_votes += taxid_votes
            if taxid_votes > taxid_leader_votes:
                taxid_leader = taxid
                taxid_leader_votes = taxid_votes
        majority_threshold = float(total_votes)/2
        if taxid_leader_votes > majority_threshold and taxid_leader != 'unclassified':
            return taxid_leader
    # Just in case
    return 1

def parse_names(names_dmp_path):
    names = {}
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid names')
    wc_output = subprocess.check_output(['wc', '-l', names_dmp_path])
    wc_list = wc_output.split()
    number_of_lines = int(wc_list[0])
    with open(names_dmp_path) as names_dmp:
        for line in tqdm(names_dmp, total=number_of_lines):
            taxid, name, _, classification = line.strip('\t|\n').split('\t|\t')[:4]
            taxid = int(taxid)
            # Only add scientific name entries
            scientific = classification == 'scientific name'
            if scientific:
                # line_list[1] = line_list[1].replace(' ', '_')
                names.update({taxid:name})
    return(names)

def parse_nodes(nodes_dmp_path):
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ' Processing taxid nodes')
    wc_output = subprocess.check_output(['wc', '-l', nodes_dmp_path])
    wc_list = wc_output.split()
    number_of_lines = int(wc_list[0])
    nodes_dmp = open(nodes_dmp_path)
    root_line = nodes_dmp.readline()
    nodes = {}
    nodes.update({1:{'parent':1, 'rank':'root'}})
    for line in tqdm(nodes_dmp, total=number_of_lines):
        child, parent, rank = line.split('\t|\t')[:3]
        parent, child = map(int,[parent, child])
        nodes.update({child:{'parent':parent, 'rank':rank}})
    nodes_dmp.close()
    return(nodes)

def parse_lca(lca_fpath):
    print( time.strftime("%Y-%m-%d %H:%M:%S") + ' Parsing lca taxonomy table')
    # Work out number of lines in file
    wc_output = subprocess.check_output(['wc', '-l', lca_fpath])
    wc_list = wc_output.split()
    number_of_lines = int(wc_list[0])
    number_of_proteins = {}
    lca_hits = {}
    # lca_hits[contig][rank][taxid] (running total of each thing)
    fh = open(lca_fpath)
    for line in tqdm(fh, total=number_of_lines):
        orf, name, rank, taxid = line.strip().split('\t')
        contig, orf_num = orf.rsplit('_', 1)
        taxid = int(taxid)
        if taxid != 1:
            while rank not in set(CANONICAL_RANKS):
                taxid = nodes.get(taxid, {'parent':1}).get('parent')
                rank = nodes.get(taxid, {'rank':'unclassified'}).get('rank')
        # Count number of proteins per contig
        if contig in number_of_proteins:
            number_of_proteins[contig] += 1
        else:
            number_of_proteins[contig] = 1

        # Keep running total of taxids for each contig
        if contig not in lca_hits:
            lca_hits.update({contig:{rank:{taxid:1}}})
            continue

        if rank not in lca_hits[contig]:
            lca_hits[contig].update({rank:{taxid:1}})
            continue

        if taxid not in lca_hits[contig][rank]:
            lca_hits[contig][rank].update({taxid:1})
        else:
            lca_hits[contig][rank][taxid] += 1
    fh.close()
    return(lca_hits)

def rank_taxids(ctg_lcas):
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ' Ranking taxids')
    n_contigs = len(ctg_lcas)
    top_taxids = {}
    for contig in tqdm(ctg_lcas, total=n_contigs):
        acceptedTaxid = None
        for rank in CANONICAL_RANKS:
            if acceptedTaxid is not None:
                break
            # Order in descending order of votes
            if rank in ctg_lcas[contig]:
                ordered_taxids = sorted(ctg_lcas[contig][rank], key=lambda tid:ctg_lcas[contig][rank][tid], reverse=True)
                #sys.exit()
                for taxid in ordered_taxids:
                    if isConsistentWithOtherOrfs(taxid, rank, ctg_lcas[contig], nodes):
                        acceptedTaxid = taxid
                        break

        # If acceptedTaxid is still None at this point, there was some kind of
        # draw, so we need to find the lowest taxonomic level where there is a
        # majority
        if acceptedTaxid is None:
            acceptedTaxid = lowest_majority(ctg_lcas[contig], nodes)

        top_taxids[contig] = acceptedTaxid
    return top_taxids

def write_assigned(assigned_contigs, outfpath):
    """ Writes taxa assigned contigs to provided output tab-delimited file path
    """
    lines = 'contig\ttaxid\n'
    for contig,taxid in assigned_contigs.items():
        lines += '{}\t{}\n'.format(contig,taxid)
    fh = open(outfpath, 'w')
    fh.write(lines)
    fh.close()
    return outfpath

def resolve_taxon_paths(ctg2taxid):
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ' Resolving taxon paths')
    contig_paths = {}
    # {contig:{rank1:name1,rank2,name2},contig2:{rank1:name1,rank2:name2,...},...}
    n_contigs = len(ctg2taxid)
    for contig in tqdm(ctg2taxid, total=n_contigs):
        taxid = ctg2taxid[contig]
        if taxid == 1:
            contig_paths.update({contig:{'root':names[taxid]}})
        while taxid != 1:
            current_rank = nodes[taxid]['rank']
            if current_rank not in set(CANONICAL_RANKS):
                taxid = nodes[taxid]['parent']
                continue

            name = names[taxid]

            if contig not in contig_paths:
                contig_paths.update({contig:{current_rank:name}})
            else:
                contig_paths[contig][current_rank] = name
            taxid = nodes[taxid]['parent']

        for rank in CANONICAL_RANKS:
            if rank not in contig_paths[contig]:
                contig_paths[contig][rank] = 'unclassified'

    for contig in contig_paths:
        contig_paths[contig].pop('root')
        contig_paths[contig]['taxid'] = ctg2taxid[contig]

    return(contig_paths)

def write_taxa(ranked_ctgs, contig_table_fpath, outfpath):
    print(time.strftime("%Y-%m-%d %H:%M:%S") + ' Writing table')
    outfile = open(outfpath, 'w')
    fh = open(contig_table_fpath)
    header = fh.readline().rstrip('\n')
    header += '\t'.join([
        '','kingdom','phylum','class','order','family','genus','species','taxid'
        ]) + '\n'
    outfile.write(header)

    for line in fh:
        original_line = line.rstrip('\n')
        line_list = original_line.split('\t')
        ctg = line_list[0]
        if ctg not in ranked_ctgs:
            # In this case we fill up the record with 'unclassified'
            # - probably this results from the contig having no blast hits
            ranked_ctgs[ctg] = {}
            for rank in CANONICAL_RANKS:
                ranked_ctgs[ctg][rank] = 'unclassified'
                ranked_ctgs[ctg]['taxid'] = 'unclassified'

        new_line = [ranked_ctgs[ctg][rank] for rank in reversed(CANONICAL_RANKS)]
        new_line.insert(0, original_line)
        new_line.append(str(ranked_ctgs[ctg]['taxid']))

        outline = '\t'.join(map(str, new_line)) + '\n'
        outfile.write(outline)
    outfile.close()



if not len(sys.argv) == 4:
    usage = 'Usage: add_contig_taxonomy.py <taxdump_dirpath> <lca_tab> <outfpath>'
    exit(usage)

# contig_tab_path = sys.argv[1]
tax_table_path = sys.argv[2]
taxdump_dir_path = sys.argv[1]
output_file_path = sys.argv[3]

# Process NCBI taxdump files
name_fpath = os.path.join(taxdump_dir_path, 'names.dmp')
nodes_fpath = os.path.join(taxdump_dir_path, 'nodes.dmp')

# Build taxid tree structure with associated canoncial ranks and names
names = parse_names(name_fpath)
nodes = parse_nodes(nodes_fpath)

# retrieve lca taxids for each contig
classifications = parse_lca(tax_table_path)

# Vote for majority lca taxid from contig lca taxids
assigned_contigs = rank_taxids(classifications)
outfpath = write_assigned(assigned_contigs, output_file_path)
print(f'written: {outfpath}')
# Add all corresponding canonical ranks from voted taxid
# ranked_contigs = resolve_taxon_paths(assigned_contigs)

# Removing root for writing columns to table
# CANONICAL_RANKS.remove('root')

# Add assigned canonical ranks and voted taxid to contig table
# write_taxa(ranked_contigs, contig_tab_path, output_file_path)

# print('written: {}'.format(os.path.abspath(output_file_path)))
