#!/usr/bin/env python
"""
Name: Lowest Common Ancestor Script (LCA)
Function: Finds LCA of all orfs from BLAST query
Usage: For usage information, type "python lca.py -help" in the command line prompt
Input Parameters: nodes.dmp, names.dmp, prot.accession2taxid, (blast output file)
Output: (blast.lca) LCA between filtered orfs (>90% top bit score) in BLAST.tab
Author: Evan R. Rees
Required Extensions: lca_functions.so | lca_functions.c (lca_functions.pyx compiled using Cython)
"""

import time, re, argparse, os
from functools import reduce
from itertools import chain
from tqdm import tqdm
from sys import argv
from sys import exit

try:
    import lca_functions
except ImportError as failed_import:
    print("\nlca.py needs access to the cython compiled file: lca_functions.c or lca_functions.so.\n\
To compile the lca_functions.c file, navigate to the directory containing lca_functions.pyx\n\
Enter the following into the command line prompt:\n\n\
cmd line:\tpython setup_lca_functions.py build_ext --inplace\n")
    exit()

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range between 0 and 1. Must convert the percentage to a decimal" % (x,))
    return (x)

def readable_dir(prospective_dir):
  if not os.path.isdir(prospective_dir):
    raise Exception("readable_dir:{0} is not a valid path".format(prospective_dir))
  if os.access(prospective_dir, os.R_OK):
    return prospective_dir
  else:
    raise Exception("readable_dir:{0} is not a readable dir".format(prospective_dir))

parser = argparse.ArgumentParser(description="Script to find the Lowest Common Ancestor for each ORF from a BLAST output table",
epilog='''The LCA analysis output will be directed to run_taxonomy.py\n\n\

NOTE:\nLCA analysis will produce best results when database files are up to date.\n\
Database files can be automatically updated before performing LCA analysis by specifying:\n\
\"lca.py [-f] [-v] [-fail_info] database_directory <path_to_database_directory> -update BLAST output\"\n\n\

Up to date versions of nodes.dmp and names.dmp may be found at:\nftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz\n\n\
prot.accession2taxid is updated weekly and may be found at:\nftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz''',
 formatter_class=argparse.RawDescriptionHelpFormatter)

subparsers = parser.add_subparsers(title='Method for Accession of Database Files',\
help="Additional help may be provided by specifying (database_files|database_directory) -h")

parser_databasefiles = subparsers.add_parser("database_files",\
help="Accesses database files provided by individually listing nodes.dmp, names.dmp and accession2taxid")
parser_databasefiles.add_argument("nodes.dmp", help="Path to nodes.dmp file")
parser_databasefiles.add_argument("names.dmp", help="Path to names.dmp file")
parser_databasefiles.add_argument("accession2taxid", help="Path to prot.accession2taxid file")
parser_databasefiles.set_defaults(parser_databasedirectory=False, parser_databasefiles=True)

parser_databasedirectory = subparsers.add_parser("database_directory", help="Accesses database files from directory provided")
parser_databasedirectory.add_argument('path_to_database_directory', action='store',\
help='Path to directory with database files.', type=readable_dir)
parser_databasedirectory.add_argument("-update", required=False, action='store_true',\
help='Updates nodes.dmp, names.dmp and accession2taxid files within the specified directory')
parser_databasedirectory.set_defaults(parser_databasefiles=False, parser_databasedirectory=True)

parser.add_argument("-v", "--verbose", required=False, help="Indicate verbose progress reporting", action='store_true')
parser.add_argument("blast", metavar='BLAST output', help="Path to BLAST output file.")
parser.add_argument("-f", metavar='bitscore filter', required=False, default=0.9, type=restricted_float,\
help="Filter to parse percentage of top BLAST hits based on bitscore.")
parser.add_argument("-fail_info", required=False, help="Writes out files with failure taxid/orf information", action='store_true')

args = vars(parser.parse_args())

if args['parser_databasefiles']:
    nodes_path = args['nodes.dmp']
    names_path = args['names.dmp']
    accession2taxid_file = args['accession2taxid']
elif args['parser_databasedirectory']:
    taxdump_dir_path = args['path_to_database_directory'].rstrip('/')
    accession2taxid_path = taxdump_dir_path + '/prot.accession2taxid'
    names_path = taxdump_dir_path + '/names.dmp'
    nodes_path = taxdump_dir_path + '/nodes.dmp'
    accession2taxid_file = accession2taxid_path

blast_file = args['blast']
bitscore_filter = args['f']
verbose = args['verbose']
failure_tracking = args['fail_info']
taxdump_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
accession2taxid_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"

if args['update']:
    os.system('wget -v %s -O %s/taxdump.tar.gz' % (taxdump_url, taxdump_dir_path))
    os.system('tar -xvzf %s/taxdump.tar.gz -C %s names.dmp nodes.dmp' % (taxdump_dir_path, taxdump_dir_path))
    os.system('rm %s/taxdump.tar.gz' % taxdump_dir_path)
    os.system('wget -v %s -O %s.gz' % (accession2taxid_url, accession2taxid_path))
    print("Gunzipping prot.accession2taxid gzipped file\nThis may take some time...")
    os.system('gunzip -9vNf %s.gz > %s' % (accession2taxid_path, accession2taxid_path))
    print("\nnodes.dmp, names.dmp and prot.accession2taxid files updated in %s\nResuming..." % taxdump_dir_path)

start_time = time.time()
output_filename = str(blast_file.split("/")[-1].rstrip(".tab"))

if verbose:
    # Parse nodes file
    parents = dict()
    children = dict()
    taxids = dict()

    total_nodes = lca_functions.wc(nodes_path)
    print('Beginning Parsing Node Dump File')

    with open(nodes_path) as nodes:
        for i,line in tqdm(enumerate(nodes), total=total_nodes, desc='Building Tax ID Paths', dynamic_ncols=True):
            if i > 0: # We skip the first line because it says that root is its own child!
                no_whitespace_string = line.replace(' ', '')
                no_whitespace_string = no_whitespace_string.replace('\t', '')
                line_list = no_whitespace_string.split('|')
                child = int(line_list[0])
                parent = int(line_list[1])

                # Add information to datastructures
                taxids[child] = 1

                parents[child] = parent

                if parent in children:
                    children[parent].add(child)
                else:
                    children[parent] = set([child])

    print('%d Tax ID Paths Created'% total_nodes)


    #data structures for tree traversal w/ distance from root and first occurrence attributes
    tour = list()
    first_node = (0, 1, 'b')
    tour.append(first_node)
    direction = 'forward'
    dist = 0
    level = list()
    level.append(dist)

    #Traversing tree by eulerian tour
    with tqdm(total=total_nodes, desc='Constructing Tree', leave=True, dynamic_ncols=True) as treebar:
        while taxids:
            if direction == 'forward':
                # Looking for a child of tour[-1][1], will take first arbitrary one
                parent = tour[-1][1]
                if parent in children: # i.e. if parent has children
                    child = children[parent].pop()
                    new_node = (parent, child, 'b')
                    tour.append(new_node)
                    dist += 1
                    level.append(dist)
                    # Delete the child taxid from the taxids dictionary
                    taxids.pop(child, None)
                    treebar.update(1)

                    # If the set is now empty, we need to delete the parent key in children
                    if not children[parent]:
                        children.pop(parent, None)
                else:
                    direction = 'reverse'
                    # Do nothing else
            elif direction == 'reverse':
                child = tour[-1][1]
                if child in parents: # i.e. child has parents left
                    parent = parents[child]
                    new_node = (child, parent, 'l')
                    tour.append(new_node)
                    dist -= 1
                    level.append(dist)
                    # Delete the child from the parents dictionary
                    parents.pop(child, None)

                    # If parent still has children, reverse direction
                    if parent in children:
                        direction = 'forward'
                else:
                    direction = 'forward'

    elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))
    print("Elapsed time is {} (HH:MM:SS)\nFinished Traversing Tree".format(elapsed_time))

    occurrence = dict()
    tour_length=len(tour)
    #Building first occurrence dictionary
    for index, node in tqdm(enumerate(tour), leave=True, total=tour_length, desc='Constructing First Occurrence Dictionary', dynamic_ncols=True):
        child = node[1]
        if child not in occurrence:
            occurrence[child]=int(index)
        else:
            pass


    elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))
    print("Elapsed time is {} (HH:MM:SS)\nFinished building first occurrence dictionary".format(elapsed_time))


    #pp.pprint(tour)
    print("\t\tLengths:\ntour:\t\t%d\nlevel:\t\t%d\noccurrence:\t%d" % (len(tour),len(level),len(occurrence)))

    """
    Next Module: Constructs Sparse Table (Preprocessing for RMQ and LCA)
    Uses level array constructed from eulerian tour
    """

    print("Building Sparse Table")
    sparse_table = lca_functions.verbosePreprocess(level)


    elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))
    print("Elapsed time is %s (HH:MM:SS)\nSparse table created." % elapsed_time)


    """
    Constructs dictionary of {'ORF1': [accession number/accession number.version, ...], ...}
    taking accession numbers for (default 90% of) topbitscore and above
    Operating under the assumption bitscore is descending from highest to lowest for each gene
    """


    print("Parsing BLAST Output File")
    num_hits = lca_functions.wc(blast_file)
    blast_orfs = lca_functions.verboseExtract_blast(blast_file, bitscore_filter)
    num_extracted_accession_numbers = len(list(chain.from_iterable(blast_orfs.values())))
    print("Finished Building Accession Number Dictionary from BLAST file\nAfter Filtering: Extracted %d out of %d hits" % (num_extracted_accession_numbers, num_hits))

    elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))
    print("Elapsed time is {} (HH:MM:SS)\nAccession Number Set Ready for Conversion".format(elapsed_time))


    """
    Next Module: Reads in Genbank accession2taxid_file
    Converts accession numbers from blast output to tax ids in preparation for LCA algorithm
    """

    accession2taxid_dict = lca_functions.verboseProcess_accession2taxid_file(accession2taxid_file, blast_orfs)


    blast_taxids = lca_functions.verboseConvert_accession2taxid(accession2taxid_dict, blast_orfs)


    elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))
    print("Elapsed time is {} (HH:MM:SS)\nConversion from accession number to tax ID completed".format(elapsed_time))


    """
    Next Module: Performs LCA algorithm on taxids from converted BLAST accession numbers
    Uses reduce and finds LCA from list of taxids from blast output
    """


    print("Beginning LCA")
    failed_orfs = list()
    failed_taxids = list()
    lca_dict = dict()
    total_rmq = len(blast_taxids)
    for orf, taxset in tqdm(blast_taxids.iteritems(), desc="Performing RMQ", leave=True, total=total_rmq):
        lca = False
        while lca == False:
            if len(taxset) >= 2:
                try:
                    lca = reduce(lambda x,y: lca_functions.RangeMinQuery(node1=x, node2=y, tree=tour, sparse_table=sparse_table, level_array=level, first_occurrence_index=occurrence), taxset)
                    lca_dict[orf] = {'lca' : int(lca)}
                except KeyError as failed_taxid:
                    #Will raise KeyError if taxid not in tree
                    failed_taxid = str(failed_taxid)
                    failed_taxid = int(failed_taxid)
                    failure_info = (orf, failed_taxid)
                    failed_taxids.append(failure_info)
                    blast_taxids[orf].remove(failed_taxid)
                except ValueError as failed_orf:
                    failure_info = (orf, taxset)
                    failed_orfs.append(failure_info)
                    lca = True
            elif len(taxset) == 1:
                failed_taxid = int(taxset.pop())
                failure_info = (orf, failed_taxid)
                failed_taxids.append(failure_info)
                lca = failed_taxid
                lca_dict[orf] = {'lca' : lca}
            else:
                failed_taxid = taxset
                failure_info = (orf, failed_taxid)
                failed_taxids.append(failure_info)
                lca_dict[orf] = {'lca' : 1}
                #default lca to root
                lca = True

    print("Range Minimum Query Called %d Times" % lca_functions.RangeMinQuery.count())

    """
     and adds to failed_taxids list, then performs RMQ again
    When through tax ID list, it will test if LCA is given.
    (LCA Will always have at least a value of 1; breaking while loop)
    """

    if failure_tracking:
        if failed_taxids:
            print('Did Not Perform Range Minimum Query on %d Tax IDs' % len(failed_taxids))
            print("Tax IDs of Which No RMQ was Performed May be Found in %s_failed_taxids.lca" % output_filename)
            print("Either Accession Numbers Could Not be Found for Conversion to Associated Tax IDs\n\t\t\tOR\nOnly 1 Tax ID was Located, Thus Not Performing an LCA Query")
            print("Beginning Search for Accession Numbers Associated with LCA Tax IDs")
            with open(output_filename + "_failed_taxids.lca","w") as failed_ids_outfile:
                failed_ids_outfile.write('List of failed Taxids for LCA (orf\\\"t\"Tax ID)\n')
                for orfs, ids in failed_taxids:
                    failed_ids_outfile.write('%s\t%s\n' % (orfs,ids))
        else:
            pass

        if failed_orfs:
            print("%d ORFs' Taxsets Failed" % len(failed_orfs))
            print("Associated Taxsets May be Found in %s_failed_orfs.lca" % output_filename)
            with open(output_filename + "_failed_orfs.lca","w") as failed_orfs_outfile:
                for orfs, taxset in failed_orfs:
                    failed_orfs_outfile.write('%s\t%s\n' % (orfs, taxset))
        else:
            pass

    elapsed_time = time.strftime('%H:%M:%S', time.gmtime(round((time.time() - start_time),2)))
    print("Elapsed time is %s (HH:MM:SS)\nFinished Performing LCA" % elapsed_time)
    print("Preparing %s.lca" % output_filename)

    reference_taxids = dict() # {taxid:{'name':'scientific name','parent':'parent taxid', 'rank':'given rank'}, taxid2:{...},...}
    total_names = lca_functions.wc(names_path)
    with open(names_path,"r") as names_dmp:
        for line in tqdm(names_dmp,  total=total_names, desc="Adding Names to Reference Tax ID Dictionary", leave=True, dynamic_ncols=True):
            line_list = line.rstrip('\n').split('|')
            for i,value in enumerate(line_list):
                line_list[i] = value.strip()
                if line_list[3] == 'scientific name':
                    reference_taxids[int(line_list[0])] = {'name':line_list[1]}


    with open(nodes_path, "r") as nodes_dmp:
        for line in tqdm(nodes_dmp, total=total_nodes, desc="Adding Ranks to Reference Tax ID Dictionary", leave=True, dynamic_ncols=True):
            line_list = line.rstrip('\n').split('|')
            # Remove trailing and leading spaces
            for i,value in enumerate(line_list):
                line_list[i] = value.strip()
            reference_taxids[int(line_list[0])]['rank'] = line_list[2]


    total_lca_keys = len(lca_dict)
    for orf in tqdm(lca_dict.iterkeys(), total=total_lca_keys, desc="Adding LCA Rank and Name to ORFs", dynamic_ncols=True):
        lca_dict[orf]['rank'] = dict()
        lca_dict[orf]['name'] = dict()
        if lca_dict[orf]['lca'] in reference_taxids.viewkeys():
            rank = reference_taxids[lca_dict[orf]['lca']]['rank']
            name = reference_taxids[lca_dict[orf]['lca']]['name']
            lca_dict[orf]['rank'] = rank
            lca_dict[orf]['name'] = name
        else:
            lca_dict[orf]['rank'] = 'no rank'
            lca_dict[orf]['name'] = 'root'



    print("Writing out %s.lca" % output_filename)
    with open(output_filename + ".lca", "w") as lca_outfile:
        for orf in lca_dict.keys():
            rank = lca_dict[orf]['rank']
            name = lca_dict[orf]['name']
            #lca = lca_dict[orf]['lca']
            lca_outfile.write('%s\t%s\t%s\n' % (orf, name, rank))
else:
    # Parse nodes file
    parents = dict()
    children = dict()
    taxids = dict()
    with open(nodes_path) as nodes:
        for i,line in enumerate(nodes):
            if i > 0: # We skip the first line because it says that root is its own child!
                no_whitespace_string = line.replace(' ', '')
                no_whitespace_string = no_whitespace_string.replace('\t', '')
                line_list = no_whitespace_string.split('|')
                child = int(line_list[0])
                parent = int(line_list[1])

                # Add information to datastructures
                taxids[child] = 1

                parents[child] = parent

                if parent in children:
                    children[parent].add(child)
                else:
                    children[parent] = set([child])


    #data structures for tree traversal w/ distance from root and first occurrence attributes
    tour = list()
    first_node = (0, 1, 'b')
    tour.append(first_node)
    direction = 'forward'
    dist = 0
    level = list()
    level.append(dist)

    #Traversing tree by eulerian tour
    while taxids:
        if direction == 'forward':
            # Looking for a child of tour[-1][1], will take first arbitrary one
            parent = tour[-1][1]
            if parent in children: # i.e. if parent has children
                child = children[parent].pop()
                new_node = (parent, child, 'b')
                tour.append(new_node)
                dist += 1
                level.append(dist)
                # Delete the child taxid from the taxids dictionary
                taxids.pop(child, None)

                # If the set is now empty, we need to delete the parent key in children
                if not children[parent]:
                    children.pop(parent, None)
            else:
                direction = 'reverse'
                # Do nothing else
        elif direction == 'reverse':
            child = tour[-1][1]
            if child in parents: # i.e. child has parents left
                parent = parents[child]
                new_node = (child, parent, 'l')
                tour.append(new_node)
                dist -= 1
                level.append(dist)
                # Delete the child from the parents dictionary
                parents.pop(child, None)

                # If parent still has children, reverse direction
                if parent in children:
                    direction = 'forward'
            else:
                direction = 'forward'

    occurrence = dict()
    tour_length=len(tour)
    #Building first occurrence dictionary
    for index, node in enumerate(tour):
        child = node[1]
        if child not in occurrence:
            occurrence[child]=int(index)
        else:
            pass


    """
    Next Module: Constructs Sparse Table (Preprocessing for RMQ and LCA)
    Uses level array constructed from eulerian tour
    """

    sparse_table = lca_functions.Preprocess(level)

    """
    Constructs dictionary of {'ORF1': [accession number/accession number.version, ...], ...}
    taking accession numbers for (default 90% of) topbitscore and above
    Operating under the assumption bitscore is descending from highest to lowest for each gene
    """

    blast_orfs = lca_functions.Extract_blast(blast_file, bitscore_filter)

    """
    Next Module: Reads in Genbank accession2taxid_file
    Converts accession numbers from blast output to tax ids in preparation for LCA algorithm
    """

    accession2taxid_dict = lca_functions.Process_accession2taxid_file(accession2taxid_file, blast_orfs)


    blast_taxids = lca_functions.Convert_accession2taxid(accession2taxid_dict, blast_orfs)

    """
    Next Module: Performs LCA algorithm on taxids from converted BLAST accession numbers
    Uses reduce and finds LCA from list of taxids from blast output
    """
    failed_orfs = list()
    failed_taxids = list()
    lca_dict = dict()
    for orf, taxset in blast_taxids.iteritems():
        lca = False
        while lca == False:
            if len(taxset) >= 2:
                try:
                    lca = reduce(lambda x,y: lca_functions.RangeMinQuery(node1=x, node2=y, tree=tour, sparse_table=sparse_table, level_array=level, first_occurrence_index=occurrence), taxset)
                    lca_dict[orf] = {'lca' : int(lca)}
                except KeyError as failed_taxid:
                    #Will raise KeyError if taxid not in tree
                    failed_taxid = str(failed_taxid)
                    failed_taxid = int(failed_taxid)
                    failure_info = (orf, failed_taxid)
                    failed_taxids.append(failure_info)
                    blast_taxids[orf].remove(failed_taxid)
                except ValueError as failed_orf:
                    failure_info = (orf, taxset)
                    failed_orfs.append(failure_info)
                    lca = True
            elif len(taxset) == 1:
                failed_taxid = int(taxset.pop())
                failure_info = (orf, failed_taxid)
                failed_taxids.append(failure_info)
                lca = failed_taxid
                lca_dict[orf] = {'lca' : lca}
            else:
                failed_taxid = taxset
                failure_info = (orf, failed_taxid)
                failed_taxids.append(failure_info)
                lca_dict[orf] = {'lca' : 1}
                #default lca to root
                lca = True


    """
     and adds to failed_taxids list, then performs RMQ again
    When through tax ID list, it will test if LCA is given.
    (LCA Will always have at least a value of 1; breaking while loop)
    """
    if failure_tracking:
        if failed_taxids:
            with open(output_filename + "_failed_taxids.lca", "w") as failed_ids_outfile:
                failed_ids_outfile.write('List of failed Taxids for LCA (orf\\\"t\"Tax ID)\n')
                for orfs, ids in failed_taxids:
                    failed_ids_outfile.write('%s\t%s\n' % (orfs,ids))
        else:
            pass

        if failed_orfs:
            with open(output_filename + "_failed_orfs.lca", "w") as failed_orfs_outfile:
                for orfs, taxset in failed_orfs:
                    failed_orfs_outfile.write('%s\t%s\n' % (orfs, taxset))
        else:
            pass

    reference_taxids = dict() # {taxid:{'name':'scientific name','parent':'parent taxid', 'rank':'given rank'}, taxid2:{...},...}
    total_names = lca_functions.wc(names_path)
    with open(names_path,"r") as names_dmp:
        for line in names_dmp:
            line_list = line.rstrip('\n').split('|')
            for i,value in enumerate(line_list):
                line_list[i] = value.strip()
                if line_list[3] == 'scientific name':
                    reference_taxids[int(line_list[0])] = {'name':line_list[1]}

    with open(nodes_path, "r") as nodes_dmp:
        for line in nodes_dmp:
            line_list = line.rstrip('\n').split('|')
            # Remove trailing and leading spaces
            for i,value in enumerate(line_list):
                line_list[i] = value.strip()
            reference_taxids[int(line_list[0])]['rank'] = line_list[2]

    for orf in lca_dict.iterkeys():
        lca_dict[orf]['rank'] = dict()
        lca_dict[orf]['name'] = dict()
        if lca_dict[orf]['lca'] in reference_taxids.viewkeys():
            rank = reference_taxids[lca_dict[orf]['lca']]['rank']
            name = reference_taxids[lca_dict[orf]['lca']]['name']
            lca_dict[orf]['rank'] = rank
            lca_dict[orf]['name'] = name
        else:
            lca_dict[orf]['rank'] = 'no rank'
            lca_dict[orf]['name'] = 'root'

    with open(output_filename + ".lca", "w") as lca_outfile:
        for orf in lca_dict.keys():
            rank = lca_dict[orf]['rank']
            name = lca_dict[orf]['name']
            #lca = lca_dict[orf]['lca']
            lca_outfile.write('%s\t%s\t%s\n' % (orf, name, rank))
