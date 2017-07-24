"""
Name: Lowest Common Ancestor Script (LCA)
Function: Cython compiled functions for LCA analysis of contigs from BLAST query output file
Usage: lca.py {taxid reference table} {blast.tab output} {genbank prot.accession_num2taxid}
Author: Evan R. Rees
Date: {\} TBD (Launch date)
"""
import re
import numpy as np
from subprocess import check_output
from tqdm import tqdm
from itertools import chain


def wc(filename):
    "Returns number of lines in file"
    return int(check_output(["wc", "-l", filename]).split()[0])

class countcalls(object):
   "Decorator that keeps track of the number of times a function is called."

   __instances = {}

   def __init__(self, f):
      self.__f = f
      self.__numcalls = 0
      countcalls.__instances[f] = self

   def __call__(self, *args, **kwargs):
      self.__numcalls += 1
      return self.__f(*args, **kwargs)

   def count(self):
      "Return the number of times the function f was called."
      return countcalls.__instances[self.__f].__numcalls

   @staticmethod
   def counts():
      "Return a dict of {function: # of calls} for all registered functions."
      return dict([(f.__name__, countcalls.__instances[f].__numcalls) for f in countcalls.__instances])

def Preprocess(level_array):
    "Returns a sparse table with the level array associated with its respective eulerian tour from tree construction"
    n = int(len(level_array))
    num_columns = int(np.floor(np.log2(n))+1)
    sparse_table = np.empty((n,num_columns))
    #height of sparse table is from 0 to n
    #width of sparse table is from 0 to logn
    sparse_table[:,0] = [index for index,values in enumerate(level_array)]
    #Instantiates sparse table first column with indices of level array
    for col in range(1, num_columns):
        #goes from col 1 to width of array because 0th column instantiated from above
        for row in range(0, n):
            #goes from row 0 to last row (n) then will increment column +1
            if 2**col <= n:
                #assesses whether the range of indices in column will exceed the number of rows in sparse table
                if row+(2**col)-1 < n:
                    #assesses whether range of indices in column will exceed the number of rows in sparse table
                    if level_array[int(sparse_table[row, (col-1)])] < level_array[int(sparse_table[(row + 2**(col-1)), (col-1)])]:
                        #if array[index at sparse table location] < array[index at sparse table location]
                        sparse_table[row, col] = sparse_table[row, col-1]
                        #assigns index from specified sparse table location to current position in sparse table
                    else:
                        sparse_table[row, col] = sparse_table[row + (2**(col-1)), col-1]
                        #if array[index at sparse table location] >= array[index at sparse table location]
                        #assigns index from specified sparse table location to current position in sparse table
                else:
                    sparse_table[row, col] = False
                    #places zeros in positions not utilized in sparse table
    return(sparse_table)


def verbosePreprocess(level_array):
    "Returns a sparse table with the level array associated with its respective eulerian tour from tree construction"
    n = int(len(level_array))
    num_columns = int(np.floor(np.log2(n))+1)
    sparse_table = np.empty((n,num_columns))
    #height of sparse table is from 0 to n
    #width of sparse table is from 0 to logn
    sparse_table[:,0] = [index for index,values in enumerate(level_array)]
    #Instantiates sparse table first column with indices of level array
    for col in tqdm(range(1, num_columns), total=num_columns-1, desc="Cols Created", leave=True, dynamic_ncols=True ):
        #goes from col 1 to width of array because 0th column instantiated from above
        for row in tqdm(range(0, n), total=n, desc="Rows Created", leave=True, dynamic_ncols=True):
            #goes from row 0 to last row (n) then will increment column +1
            if 2**col <= n:
                #assesses whether the range of indices in column will exceed the number of rows in sparse table
                if row+(2**col)-1 < n:
                    #assesses whether range of indices in column will exceed the number of rows in sparse table
                    if level_array[int(sparse_table[row, (col-1)])] < level_array[int(sparse_table[(row + 2**(col-1)), (col-1)])]:
                        #if array[index at sparse table location] < array[index at sparse table location]
                        sparse_table[row, col] = sparse_table[row, col-1]
                        #assigns index from specified sparse table location to current position in sparse table
                    else:
                        sparse_table[row, col] = sparse_table[row + (2**(col-1)), col-1]
                        #if array[index at sparse table location] >= array[index at sparse table location]
                        #assigns index from specified sparse table location to current position in sparse table
                else:
                    sparse_table[row, col] = False
                    #places zeros in positions not utilized in sparse table
    return(sparse_table)


def Extract_blast(blast_file, bitscore_filter=0.9):
    "Returns a dictionary of accession numbers extracted from blast output file\ndefault bitscore filter takes orfs from >90% of top bitscore"
    with open(blast_file,"r") as blast_file:
        blast_dict = dict()
        temp_orf_list = list()
        for line in blast_file:
            match = re.search(r'(NODE_\d*_length_\d*_cov_\d*\.\d*_ID_\d*_\d*)\s(\S+)\s\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+\.\d)',line)
            if match:
                orf = match.group(1)
                if orf in temp_orf_list:
                    bitscore = float(match.group(3))
                    if bitscore >= bitscore_filter*topbitscore:
                        accession_num = match.group(2)
                        blast_dict[orf].append(accession_num)
                else:
                    accession_num = match.group(2)
                    bitscore = float(match.group(3))
                    blast_dict[orf] = [accession_num]
                    topbitscore = bitscore
                    temp_orf_list = []
                    temp_orf_list.append(orf)
    return(blast_dict)

def verboseExtract_blast(blast_file, bitscore_filter=0.9):
    "Returns a dictionary of accession numbers extracted from blast output file\ndefault bitscore filter takes orfs from >90% of top bitscore"
    total_hits = wc(blast_file)
    with open(blast_file,"r") as blast_file:
        blast_dict = dict()
        temp_orf_list = list()
        for line in tqdm(blast_file, desc='Filtering Accession Numbers', leave=True, total=total_hits, dynamic_ncols=True):
            match = re.search(r'(NODE_\d*_length_\d*_cov_\d*\.\d*_ID_\d*_\d*)\s(\S+)\s\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+\.\d)',line)
            if match:
                orf = match.group(1)
                if orf in temp_orf_list:
                    bitscore = float(match.group(3))
                    if bitscore >= bitscore_filter*topbitscore:
                        accession_num = match.group(2)
                        blast_dict[orf].append(accession_num)
                else:
                    accession_num = match.group(2)
                    bitscore = float(match.group(3))
                    blast_dict[orf] = [accession_num]
                    topbitscore = bitscore
                    temp_orf_list = []
                    temp_orf_list.append(orf)
    return(blast_dict)


def benchmarking_extract_blast(blast_file, bitscore_filter=0.9):
    "Returns a dictionary of accession numbers extracted from blast output file\ndefault bitscore filter takes orfs from >90% of top bitscore"
    total_hits = wc(blast_file)
    with open(blast_file,"r") as blast_file:
        blast_dict = dict()
        temp_list = list()
        for line in tqdm(blast_file, desc='Filtering Accession Numbers', leave=True, total=total_hits, dynamic_ncols=True):
            match = re.search(r'(\S+)\s\S+ref\|(\S+)\|\s\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+\.\d)',line)
            if match:
                orf = match.group(1)
                if orf in temp_list:
                    bitscore = float(match.group(3))
                    if bitscore >= bitscore_filter*topbitscore:
                        accession_num = match.group(2)
                        blast_dict[orf].append(accession_num)
                else:
                    accession_num = match.group(2)
                    bitscore = float(match.group(3))
                    blast_dict[orf]=[accession_num]
                    topbitscore = bitscore
                    temp_list = []
                    temp_list.append(orf)
    return(blast_dict)

def Convert_accession2taxid(accession2taxid_dict, blast_dict):
    "Returns dictionary of orfs with set of their associated tax ids converted from accesssion numbers and accession.version numbers"
    blast_taxids = dict()
    for orf in blast_dict.iterkeys():
        blast_taxids[orf]=set()
        for accession_number in blast_dict[orf].__iter__():
            if accession_number in accession2taxid_dict.viewkeys():
                blast_taxids[orf].add(int(accession2taxid_dict.get(accession_number)))
    return(blast_taxids)

def verboseConvert_accession2taxid(accession2taxid_dict, blast_dict):
    "Returns dictionary of orfs with set of their associated tax ids converted from accesssion numbers and accession.version numbers"
    blast_taxids = dict()
    for orf in tqdm(blast_dict.iterkeys(), desc='Converting Orfs', leave=True, dynamic_ncols=True):
        blast_taxids[orf]=set()
        for accession_number in tqdm(blast_dict[orf].__iter__(), desc='Converting Orf Accession Numbers to Tax IDs', leave=True, dynamic_ncols=True):
            if accession_number in accession2taxid_dict.viewkeys():
                blast_taxids[orf].add(int(accession2taxid_dict.get(accession_number)))
    return(blast_taxids)

def Process_accession2taxid_file(accession2taxid_file, blast_dict):
    "Returns dictionary of accession numbers and accession.version with their associated tax ids"
    accession2taxid_dict = dict()
    accession_set = set(chain.from_iterable(blast_dict.values()))
    with open(accession2taxid_file,"r") as accession2taxid_file:
        for line in accession2taxid_file:
            pattern = re.search(r'(\S+)\s+(\S+)\s+(\d+)\s+\S+', line)
            if pattern:
                accession_number = pattern.group(1)
                accession_version = pattern.group(2)
                if accession_number in accession_set:
                    taxid = int(pattern.group(3))
                    accession2taxid_dict[accession_number]=taxid
                elif accession_version in accession_set:
                    taxid = int(pattern.group(3))
                    accession2taxid_dict[accession_version]=taxid
    return(accession2taxid_dict)

def verboseProcess_accession2taxid_file(accession2taxid_file, blast_dict):
    "Returns dictionary of accession numbers and accession.version with their associated tax ids"
    accession2taxid_dict = dict()
    total_lines=wc(accession2taxid_file)
    accession_set = set(chain.from_iterable(blast_dict.values()))
    with open(accession2taxid_file,"r") as accession2taxid_file:
        for line in tqdm(accession2taxid_file, desc="Building Accession Conversion Dictionary", leave=True, total=total_lines, dynamic_ncols=True):
            pattern = re.search(r'(\S+)\s+(\S+)\s+(\d+)\s+\S+', line)
            if pattern:
                accession_number = pattern.group(1)
                accession_version = pattern.group(2)
                if accession_number in accession_set:
                    taxid = int(pattern.group(3))
                    accession2taxid_dict[accession_number]=taxid
                elif accession_version in accession_set:
                    taxid = int(pattern.group(3))
                    accession2taxid_dict[accession_version]=taxid
    return(accession2taxid_dict)


#Empty data structures used for RangeMinQuery
tour=list()
sparse_table=list()
level = list()
occurrence=dict()

@countcalls
def RangeMinQuery(node1, node2, tree=tour, sparse_table=sparse_table, level_array=level, first_occurrence_index=occurrence):
    "Returns LCA by performing Range Minimum Query b/w 2 tax ids\nRequires:euler tree, sparse table, level array, occurrence dictionary"
    if node1 == node2:
        #If after performing previous RMQs the reduced input is the same tax ID as the next input tax ID, no RMQ is needed, thus the same tax ID will be returned.
        return(node1)
    elif first_occurrence_index[node1] < first_occurrence_index[node2]:
        #From tax IDs converted from filtered BLAST output accession numbers, will find the tax ID in the first occurrence dictionary and evaluate position in tour
        #Want the tax ID that occurs first to be the low. This 'if,else' statement will appropriate tax IDs to corrrect assignment for range query.
        low = first_occurrence_index[node1]
        high = first_occurrence_index[node2]
    else:
        low = first_occurrence_index[node2]
        high = first_occurrence_index[node1]
    cutoff_range = int(np.floor(np.log2(high-low+1)))
    #The cut off range determines the two partitions in the sparse table relevant for level value lookup.
    #This allows for equipartitionining of the range b/w both nodes.
    if level_array[int(sparse_table[low, cutoff_range])] <= level_array[int(sparse_table[(high-(2**cutoff_range)+1), cutoff_range])]:
        #If the level at the index in the sparse table lower range is less than or equal to the level at index in the upper range...
        return(tree[level_array.index(level_array[int( sparse_table[low, cutoff_range] ) ], low, high)][1])
        #LCA is minimum in lower range
    else:
        #If the level at the index in the sparse table lower range is greater than the level at index in the upper range...
        return(tree[level_array.index(level_array[int( sparse_table[(high-(2**cutoff_range)+1), cutoff_range] ) ], low, high)][1])
        #LCA is minimum in upper range
