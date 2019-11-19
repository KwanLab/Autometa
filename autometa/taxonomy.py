#!/usr/bin/env python3
"""
Script to retrieve taxonomic information for metagenome dataset
"""


import argparse
import gzip
import os
import subprocess
import sys
import pickle

import numpy as np
import pandas as pd

from itertools import chain
from functools import reduce
# Non-built-in Dependencies
from Bio import SeqIO
from tqdm import tqdm
# Autometa Imports
import dependencies

class TaxonUtils(object):
    """docstring for TaxonUtils."""

    BLASTP_EXT = '.blastp'

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

    def __init__(self, dbconfig, outdir=None, verbose=False, force=False, cpus=1, usepickle=True):
        self._config = dbconfig
        self._verbose = verbose
        self._force = force
        self._tqdm = not verbose
        self._usepickle = usepickle
        self._cpus = cpus
        self._names_fpath = self._config.get('ncbi','names')
        self._nodes_fpath = self._config.get('ncbi','nodes')
        self._accession2taxid_fpath = self._config.get('ncbi','accession2taxid')
        self._blastdb_fpath = self._config.get('ncbi','blastdb')
        self._sparse_pkl = os.path.join(self._config.get('common','ncbi_dir'), 'sparse.pkl')
        self._tour = os.path.join(self._config.get('common','ncbi_dir'), 'tour.pkl')
        self._level = os.path.join(self._config.get('common','ncbi_dir'), 'level.pkl')
        self._occurrence = os.path.join(self._config.get('common','ncbi_dir'), 'occurrence.pkl')
        self._outdir = os.path.realpath(outdir) if outdir else os.path.realpath(os.curdir)
        self._pkl_fname = os.path.splitext(os.path.basename(self._accession2taxid_fpath))[0]+'.pkl'
        self._accession2taxid_pkl = os.path.join(self._outdir, self._pkl_fname)
        # Must remain below file path definitions
        self._nodes = self.parse_nodes()
        self._names = self.parse_names()
        self._lca_prepared = False
        self._acc2taxids_prepared = False
        self._taxa_assigned = False

    # ========================================
    # Definitions of Taxonomy Class Properties
    # ========================================

    def pipeline():
        doc = "The pipeline property."
        def fget(self):
            return os.path.join(os.path.dirname(os.path.dirname(__file__)), 'pipeline')
        def fset(self, value):
            self._pipeline = value
        def fdel(self):
            del self._pipeline
        return locals()
    pipeline = property(**pipeline())

    def names_fpath():
        doc = "The names_fpath property."
        def fget(self):
            return self._names_fpath
        def fset(self, value):
            self._names_fpath = value
        def fdel(self):
            del self._names_fpath
        return locals()
    names_fpath = property(**names_fpath())

    def names():
        doc = "The names property."
        def fget(self):
            return self._names
        def fset(self, value):
            self._names = value
        def fdel(self):
            del self._names
        return locals()
    names = property(**names())

    def nodes_fpath():
        doc = "The nodes_fpath property."
        def fget(self):
            return self._nodes_fpath
        def fset(self, value):
            self._nodes_fpath = value
        def fdel(self):
            del self._nodes_fpath
        return locals()
    nodes_fpath = property(**nodes_fpath())

    def nodes():
        doc = "The nodes property."
        def fget(self):
            return self._nodes
        def fset(self, value):
            self._nodes = value
        def fdel(self):
            del self._nodes
        return locals()
    nodes = property(**nodes())

    def accession2taxid_fpath():
        doc = "The accession2taxid_fpath property."
        def fget(self):
            return self._accession2taxid_fpath
        def fset(self, value):
            self._accession2taxid_fpath = value
        def fdel(self):
            del self._accession2taxid_fpath
        return locals()
    accession2taxid_fpath = property(**accession2taxid_fpath())

    def acc2taxids():
        doc = "The acc2taxids property."
        def fget(self):
            return self._acc2taxids
        def fset(self, value):
            self._acc2taxids = value
        def fdel(self):
            del self._acc2taxids
        return locals()
    acc2taxids = property(**acc2taxids())

    def blastdb_fpath():
        doc = "The blastdb_fpath property."
        def fget(self):
            return self._blastdb_fpath
        def fset(self, value):
            self._blastdb_fpath = value
        def fdel(self):
            del self._blastdb_fpath
        return locals()
    blastdb_fpath = property(**blastdb_fpath())

    def blastp_fpath():
        doc = "The blastp_fpath property."
        def fget(self):
            return self._blastp_fpath
        def fset(self, value):
            self._blastp_fpath = value
        def fdel(self):
            del self._blastp_fpath
        return locals()
    blastp_fpath = property(**blastp_fpath())

    def outdir():
        doc = "The outdir property."
        def fget(self):
            return self._outdir
        def fset(self, value):
            self._outdir = value
        def fdel(self):
            del self._outdir
        return locals()
    outdir = property(**outdir())

    def accession2taxid_pkl():
        doc = "The accession2taxid_pkl property."
        def fget(self):
            return self._accession2taxid_pkl
        def fset(self, value):
            self._accession2taxid_pkl = value
        def fdel(self):
            del self._accession2taxid_pkl
        return locals()
    accession2taxid_pkl = property(**accession2taxid_pkl())

    def blastp_hits():
        doc = "The blastp_hits property."
        def fget(self):
            return self._blastp_hits
        def fset(self, value):
            self._blastp_hits = value
        def fdel(self):
            del self._blastp_hits
        return locals()
    blastp_hits = property(**blastp_hits())

    def lca_prepared():
        doc = "The lca_prepared property."
        def fget(self):
            return self._lca_prepared
        def fset(self, value):
            self._lca_prepared = bool(value)
        def fdel(self):
            del self._lca_prepared
        return locals()
    lca_prepared = property(**lca_prepared())

    def lcas():
        doc = "The lcas property."
        def fget(self):
            return self._lcas
        def fset(self, value):
            self._lcas = value
        def fdel(self):
            del self._lcas
        return locals()
    lcas = property(**lcas())

    def taxa_assigned():
        doc = "The taxa_assigned property."
        def fget(self):
            return self._taxa_assigned
        def fset(self, value):
            self._taxa_assigned = bool(value)
        def fdel(self):
            del self._taxa_assigned
        return locals()
    taxa_assigned = property(**taxa_assigned())

    # ========================================
    # Definitions of Taxonomy Class Methods
    # ========================================

    def unpickle(self, fpath):
        """ Load a pickle file. Opposite of make_pickle method that writes
        object to a pickle file (*.pkl).

        Input:
            fpath - </path/to/pickled/file>
        ReturnType: mixed
        Returns:
            obj - python object that was serialized to file via make_pickle
        """
        # NOTE: Load pickle probably should be in a common utilities class outside of taxonomy.py
        if self._verbose:
            print(f'unpickling {fpath}')
        with open(fpath, 'rb') as fh:
            obj = pickle.load(file=fh)
        if self._verbose:
            print(f'{fpath} object unpickled')
        return obj

    def make_pickle(self, obj, outfpath):
        """ Serialize a python object to a pickle file (*.pkl). Opposite of
        unpickle method that retrieves python object from file.
        Input:
            obj - python object to serialize to outfpath
            outfpath - </path/to/pickled/file>
        ReturnType: str
        Returns:
            outfpath - </path/to/pickled/file>
        """
        if self._verbose:
            print(f'pickling object to {outfpath}')
        with open(outfpath, 'wb') as fh:
            pickle.dump(obj=obj, file=fh)
        if self._verbose:
            print(f'object pickled to {outfpath}')
        return outfpath

    def name(self, taxid, rank=None):
        # If taxid is not found in names.dmp, will return None
        try:
            taxid = int(taxid)
        except ValueError as err:
            print(f'Taxid must be an integer! {taxid} type --> {type(taxid)}')
            return None
        if not rank:
            return self._names.get(taxid, 'unclassified')
        if rank not in set(self.CANONICAL_RANKS):
            print(f'{rank} not in canonical ranks!')
            return None
        ancestor_taxid = taxid
        while ancestor_taxid != 1:
            ancestor_rank = self.rank(ancestor_taxid)
            if ancestor_rank == rank:
                return self._names.get(ancestor_taxid)
            ancestor_taxid = self.parent(ancestor_taxid)

    def lineage(self, taxid):
        """ Returns the lineage of taxids encountered when traversing to root
        Input:
            taxid - taxid in nodes.dmp (int)
        ReturnType: ordered list of dicts
        Returns: [{'taxid':taxid, rank':rank,'name':'name'}, ...]
        """
        try:
            taxid = int(taxid)
        except ValueError as err:
            print(f'Taxid must be an integer! {taxid} type --> {type(taxid)}')
            return None
        lineage = []
        while taxid != 1:
            lineage.append({
            'taxid':taxid,
            'name':self.name(taxid),
            'rank':self.rank(taxid)})
            taxid = self.parent(taxid)
        return lineage

    def rank(self, taxid):
        # If taxid is not found in nodes.dmp, will return 'unclassified'
        try:
            taxid = int(taxid)
        except ValueError as err:
            print(f'Taxid must be an integer! {taxid} type --> {type(taxid)}')
            return None
        return self._nodes.get(taxid,{'rank':'unclassified'}).get('rank')

    def parent(self, taxid):
        # If taxid is not found in nodes.dmp, will return 1
        try:
            taxid = int(taxid)
        except ValueError as err:
            print(f'Taxid must be an integer! {taxid} type --> {type(taxid)}')
            return None
        return self._nodes.get(taxid,{'parent':1}).get('parent')

    def parse_names(self):
        if self._verbose:
            print(f'Processing names from {self.names_fpath}')
        names = {}
        fh = open(self.names_fpath)
        for line in tqdm(fh, disable=self._tqdm, desc='parsing names', leave=False):
            taxid, name, _, classification = line.strip('\t|\n').split('\t|\t')[:4]
            taxid = int(taxid)
            name = name.lower()
            # Only add scientific name entries
            is_scientific = classification == 'scientific name'
            if is_scientific:
                names.update({taxid:name})
        fh.close()
        if self._verbose:
            print('names loaded')
        return(names)

    def parse_nodes(self):
        if self._verbose:
            print(f'Processing nodes from {self.nodes_fpath}')
        fh = open(self.nodes_fpath)
        __ = fh.readline() # root line
        nodes = {1:{'parent':1, 'rank':'root'}}
        for line in tqdm(fh, disable=self._tqdm, desc='parsing nodes', leave=False):
            child, parent, rank = line.split('\t|\t')[:3]
            parent, child = map(int,[parent, child])
            rank = rank.lower()
            nodes.update({child:{'parent':parent,'rank':rank}})
        fh.close()
        if self._verbose:
            print('nodes loaded')
        return(nodes)

    def __prepare_tree(self):
        if self._usepickle:
            repickle = False
            # Creating a repickle toggle b/c if *any* of the fpaths do not exist need to redo
            for fp in [self._tour, self._level, self._occurrence]:
                if not os.path.exists(fp):
                    repickle = True
            if not repickle:
                self.__tour = self.unpickle(fpath=self._tour)
                self.__level = self.unpickle(fpath=self._level)
                self.__occurrence = self.unpickle(fpath=self._occurrence)
                return

        if self._verbose:
            print('Preparing tree, level, occurrence for LCA/RMQ')
        taxids = {}
        parents = {}
        children = {}
        with open(self.nodes_fpath) as fh:
            _ = fh.readline() #root
            for line in fh:
                child, parent = line.split('\t|\t')[:2]
                parent, child = map(int,[parent,child])
                taxids.update({child:1})
                parents.update({child:parent})
                if parent in children:
                    children[parent].add(child)
                else:
                    children.update({parent:set([child])})
        self.__tour = [(0,1)]
        direction = 1
        dist = 0
        self.__level = [dist]
        while taxids:
            if direction > 0:
                parent = self.__tour[-1][1]
                if parent not in children:
                    direction *= -1
                    continue
                child = children[parent].pop()
                new_node = (parent, child)
                self.__tour.append(new_node)
                dist += 1
                self.__level.append(dist)
                taxids.pop(child, None)
                # If the set is now empty, we need to delete the parent key in children
                if not children[parent]:
                    children.pop(parent, None)
                # Do nothing else
            elif direction < 0:
                child = self.__tour[-1][1]
                if child not in parents:
                    direction *= -1
                    continue
                parent = parents[child]
                new_node = (child, parent)
                self.__tour.append(new_node)
                dist -= 1
                self.__level.append(dist)
                # Delete the child from the parents dictionary
                parents.pop(child, None)
                # If parent still has children, reverse direction
                if parent in children:
                    direction *= -1
        self.__occurrence = {}
        for i,node in enumerate(self.__tour):
            child = node[1]
            if child not in self.__occurrence:
                self.__occurrence[child]=i
        if self._usepickle:
            self.make_pickle(obj=self.__tour, outfpath=self._tour)
            self.make_pickle(obj=self.__level, outfpath=self._level)
            self.make_pickle(obj=self.__occurrence, outfpath=self._occurrence)

    def __preprocess_minimums(self):
        """
        Returns:
            sparse table of minimums to be used for LCA/RangeMinimumQuery
        Uses level array associated with its respective eulerian tour from tree construction
        See method __prepare_tree
        sparse table:
            n = number of elements in levels list
            nrows - (0 - n)
            ncols - (0 - logn)
        """
        if self._usepickle and os.path.exists(self._sparse_pkl):
            self._sparse_matrix = self.unpickle(fpath=self._sparse_pkl)
            return

        if self._verbose:
            print('Constructing Sparse Table')
        nrows = len(self.__level)
        ncols = int(np.floor(np.log2(nrows))+1)
        self._sparse_matrix = np.empty((nrows,ncols))
        self._sparse_matrix[:,0] = [i for i in range(nrows)]
        for col in tqdm(range(1, ncols), disable=self._tqdm, desc='Precomputing LCAs'):
            for row in range(0, nrows):
                if 2**col > nrows:
                    continue
                if row+(2**col)-1 >= nrows:
                    self._sparse_matrix[row, col] = False
                    continue
                lower_index = self._sparse_matrix[row, (col-1)]
                upper_index = self._sparse_matrix[(row + 2**(col-1)), (col-1)]
                lower_index, upper_index = map(int, [lower_index, upper_index])
                lower_min = self.__level[lower_index]
                upper_min = self.__level[upper_index]
                if lower_min < upper_min:
                    self._sparse_matrix[row, col] = lower_index
                else:
                    self._sparse_matrix[row, col] = upper_index
        if self._usepickle:
            self.make_pickle(obj=self._sparse_matrix, outfpath=self._sparse_pkl)

    def _prepare_lca(self):
        if self._verbose:
            print('Preparing data structures for LCA')
        self.__prepare_tree()
        self.__preprocess_minimums()
        self._lca_prepared = True
        # self.__tour, self.__level, self.__occurrence, self._sparse_matrix

    def lca(self, node1, node2):
        """ Performs Range Minimum Query b/w 2 taxids
        Inputs:
            node1 - taxid (int)
            node2 - taxid (int)
        ReturnType: int
        Returns: LCA
        For requirements see method _prepare_lca()
        """
        if not self.lca_prepared:
            self._prepare_lca()
        if node1 is None and node2 is None:
            return 1
        if node1 is None:
            return node2
        if node2 is None:
            return node1
        if node1 == node2:
            return node1
        if node1 not in self.__occurrence:
            print(f'{node1} not in tree')
            return node1
        if node2 not in self.__occurrence:
            print(f'{node2} not in tree')
            return node2
        if self.__occurrence[node1] < self.__occurrence[node2]:
            low = self.__occurrence[node1]
            high = self.__occurrence[node2]
        else:
            low = self.__occurrence[node2]
            high = self.__occurrence[node1]
        # equipartition range b/w both nodes.
        cutoff_range = int(np.floor(np.log2(high-low+1)))
        lower_index = self._sparse_matrix[low, cutoff_range]
        upper_index = self._sparse_matrix[(high-(2**cutoff_range)+1), cutoff_range]
        lower_index, upper_index = map(int, [lower_index, upper_index])
        lower_range = self.__level[lower_index]
        upper_range = self.__level[upper_index]
        if lower_range <= upper_range:
            lca_range = lower_range
        else:
            lca_range = upper_range
        return self.__tour[self.__level.index(lca_range, low, high)][1]

    def get_acc2taxids(self, hits=None):
        """
        Input:
            hits - {qseqid:set(sseqids)} See return of self.parse_blastp
        ReturnType: dict
        Returns:
            self._taxids - {qseqid:{taxid, taxid, taxid, ...}, qseqid:{...}, ...}
        """
        if self._verbose:
            print(f'Converting accesssions to taxids')
        if not hits and not self._blastp_hits:
            print(f'HitsNotFoundError:\nhits:{hits}\nblastp_hits:{self._blastp_hits}')
            return None
        if not hits and self._blastp_hits:
            hits = self._blastp_hits
        if not self._acc2taxids_prepared:
            self._process_acc2taxid_db(hits=hits)
        self._taxids = {}
        old_db = False
        for qseqid,sseqids in hits.items():
            tids = {self.acc2taxids.get(sseqid) for sseqid in sseqids}
            # NOTE: Instances occur where taxid is None due to accession ID being supressed.
            # I.e. Accession number is deprecated and nr.dmnd needs to be updated.
            if None in tids:
                old_db = True
            tids.discard(None)
            self._taxids.update({qseqid:tids})
        if old_db and self._verbose:
            print(
                'WARNING: Some accessions are missing from the prot.accession2taxid'
                ' database.\nConsider updating the blast database, as these'
                ' accessions are likely now deprecated/suppressed.'
            )
        return self._taxids

    def _process_acc2taxid_db(self, hits=None):
        """ Constructs dict by searching through prot.accession2taxid.gz for
        accessions retrieved from parse_blastp
        Input:
            hits - See parse_blastp
        ReturnType: dict
        Returns: {accession:taxid, accession: taxid, ...}
        """
        if self._verbose and self._usepickle:
            print(f'Searching for {self._accession2taxid_pkl}')
        all = False
        if not hits and not self._blastp_hits:
            all = True
        if not hits and self._blastp_hits:
            all = False
            hits = self._blastp_hits
        if os.path.exists(self._accession2taxid_pkl) and self._usepickle:
            self._acc2taxids = self.unpickle(fpath=self._accession2taxid_pkl)
            self._acc2taxids_prepared = True
            return self._acc2taxids
        if self._verbose:
            print(f'Searching for {self._accession2taxid_fpath}')
        if not os.path.exists(self._accession2taxid_fpath):
            print(f'FileNotFoundError {self._accession2taxid_fpath}')
            return None
        acc2taxids = {}
        accessions = set(chain.from_iterable(self.blastp_hits.values()))
        if self._accession2taxid_fpath.endswith('.gz'):
            fh = gzip.open(self._accession2taxid_fpath)
            is_gzipped = True
        else:
            fh = open(self._accession2taxid_fpath)
            is_gzipped = False
        __ = fh.readline()
        for line in tqdm(fh, disable=self._tqdm, desc='Searching acc2taxid DB', leave=False):
            if is_gzipped:
                line = line.decode()
            acc_num, acc_ver, taxid, _ = line.split('\t')
            taxid = int(taxid)
            if acc_num in accessions:
                acc2taxids.update({acc_num:taxid})
            if acc_ver in accessions:
                acc2taxids.update({acc_ver:taxid})
        fh.close()
        self._acc2taxids = acc2taxids
        self._acc2taxids_prepared = True
        if self._usepickle:
            self.make_pickle(obj=self._acc2taxids, outfpath=self._accession2taxid_pkl)

    def acc2taxid(self, accession):
        return self._acc2taxids.get(accession, None)

    def search_blastdb(self, fasta, out=None, evalue=float('1e-5'), maxtargetseqs=200, tmpdir=os.curdir):
        """
        Will blastp search fasta file against self.blastdb_fpath from config
        Inputs:
            fasta - </path/to/protein/fasta/file>
            out - </path/to/output/file>
            evalue - (float) cutoff e-value to count hit as significant
            maxtargetseqs - (int) max number of target sequences to return
            tmpdir - </path/to/temporary/directory>
        ReturnType: str
        Returns:
            out - </path/to/output/file>
        """
        if not os.path.exists(fasta):
            raise FileNotFoundError(fasta)
        if not out:
            outdir = os.path.dirname(os.path.realpath(fasta))
            outfname = os.path.splitext(os.path.basename(fasta))[0]+self.BLASTP_EXT
            out = os.path.join(outdir, outfname)
        if os.path.exists(out) and not self._force:
            if self._verbose:
                print(f'FileAlreadyExists: {out}. To overwrite use --force')
            self._blastp_fpath = out
            return out
        if self._verbose:
            print(f'HomologySearch {fasta} against {self.blastdb_fpath}')
        cmd = map(str, [
            'diamond',
            'blastp',
            '--query',
            fasta,
            '--db',
            self.blastdb_fpath,
            '--evalue',
            evalue,
            '--max-target-seqs',
            maxtargetseqs,
            '--threads',
            self._cpus,
            '--outfmt',
            '6',
            '--out',
            out,
            '--tmpdir',
            tmpdir,
        ])
        if self._verbose:
            print(f'RunningDiamond: {" ".join(cmd)}')
        with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
            proc = subprocess.run(cmd, stdout=stdout, stderr=stderr)
        if proc.returncode:
            print(f'DiamondFailed:\nArgs:{proc.args}\nReturnCode:{proc.returncode}')
            return None
        self._blastp_fpath = out
        return out

    def parse_blastp(self, fpath=None, top_pct=0.9):
        """ Returns a dictionary of sseqids extracted from BLAST outfmt 6
        Input:
            fpath - </path/to/blastp/outfmt6/output/file>
            top_pct - bitscore filter applied to each qseqid <float in range(0,1)>
        ReturnType: dict
        Returns {qseqid:set([sseqid, sseqid,...]),...}
        """
        if not fpath:
            fpath = self.blastp_fpath
        if self._verbose:
            print(f'Parsing accessions from {fpath}')
        if not os.path.exists(fpath):
            print(f'FileNotFoundError: {fpath}')
            return None
        try:
            float(top_pct)
        except ValueError as err:
            print(f'top_pct must be a float! Input: {top_pct} Type: {type(top_pct)}')
            return None
        in_range = 0.0 < top_pct <= 1.0
        if not in_range:
            print(f'top_pct not in range(0,1)! Input: {top_pct}')
            return None
        hits = {}
        temp = set()
        with open(fpath) as fh:
            for line in tqdm(fh, disable=self._tqdm, desc='Parsing Accessions',leave=False):
                qseqid, sseqid, *__, bitscore = line.rstrip().split('\t')
                bitscore = float(bitscore)
                if qseqid not in temp:
                    hits.update({qseqid:set([sseqid])})
                    topbitscore = bitscore
                    temp = set([qseqid])
                    continue
                if bitscore >= top_pct * topbitscore:
                    hits[qseqid].add(sseqid)
        self._blastp_hits = hits
        return hits

    def get_lcas(self, taxid_sets=None):
        """ Retrieves the LCA of the taxids provided in a taxid set of arbitrary length
        Input: {qseqid:set(taxids)} For more details see return of self.get_acc2taxids method
        ReturnType: dict
        Returns: {qseqid:lca,qseqid:lca,...}
        """
        if not taxid_sets:
            taxid_sets = self._taxids
        lcas = {}
        for qseqid,taxids in taxid_sets.items():
            lca = False
            num_taxids = len(taxids)
            while not lca:
                if num_taxids >= 2:
                    lca = reduce(lambda taxid1,taxid2: self.lca(node1=taxid1,node2=taxid2), taxids)
                if num_taxids == 1:
                    lca = taxids.pop()
                # Exception handling where input for qseqid contains no taxids
                if num_taxids == 0:
                    lca = 1
            lcas.update({qseqid:lca})
        self._lcas = lcas
        return lcas

    def blast2lca(self, fasta):
        if self._verbose:
            print(f'Running BLAST to LCA for {fasta}')
        blastp_fpath = self.search_blastdb(fasta=fasta)
        hits = self.parse_blastp(fpath=blastp_fpath)
        hits = self.get_acc2taxids(hits=hits)
        return self.get_lcas(taxid_sets=hits)

    def aggregate_lcas(self, orfs_lcas=None):
        """ Aggregates ORFs' LCAS for each respective contig
        Input:
            orfs_lcas - See get_lcas method for input data structure
        ReturnType: dict
        Returns: {contig:{taxid:count, taxid2:count,...}, ...}
        """
        if self._verbose:
            print('Aggregating Contig LCAs')
        if not orfs_lcas:
            orfs_lcas = self._lcas
        contigs = {}
        for orf,taxid in orfs_lcas.items():
            contig, orf_num = orf.rsplit('_',1)
            if contig not in contigs:
                contigs.update({contig:{taxid:1}})
            elif taxid not in contigs[contig]:
                contigs[contig][taxid] = 1
            else:
                contigs[contig][taxid] += 1
        return contigs

    def assign_taxa(self, fasta, method='majority_vote'):
        """ Assigns taxa to input fasta contigs
        Inputs:
            fasta - </path/to/fasta/file>
            method - selection of how to assign taxon to a contig
        ReturnType: dict
        Returns: {contig:taxid, contig:taxid, ...}

        Methods:
        * majority_vote *
        Will determine LCAs for each ORF on the contig and find the taxid
        that is the majority leader and highest resolution by rank.
        Ranks used are the seven canonical ranks:
        superkingdom, phylum, class, order, family, genus, species.
        NOTE: Does NOT detect whether contig is from a novel taxon.

        * hierarchical *
        # TODO: See NSF Career proposal
        Raises NotImplementedError
        """
        method = method.lower()
        if self._verbose:
            print(f'Assigning taxa to {fasta} using {method} method')
        if method == 'majority_vote':
            # rank_taxids (group ORFs by contig)
            orfs_lcas = self.blast2lca(fasta)
            lcas_fpath = os.path.join(self.outdir, 'orfs_lcas.tsv')
            lcas_fpath = self.write_lcas(lcas=orfs_lcas, outfpath=lcas_fpath)
            vote_fpath = os.path.join(self.outdir, 'majority_vote.tsv')
            return self.majority_vote(infpath=lcas_fpath, outfpath=vote_fpath)

        elif method == 'hierarchical':
            contig_lcas = self.aggregate_lcas(orfs_lcas)
            raise NotImplementedError
        raise NotImplementedError

    def write_lcas(self, outfpath, lcas=None):
        # lines = 'orf\tname\trank\ttaxid\n'
        lcas = lcas if lcas else self.lcas
        lines = ''
        for orf,taxid in lcas.items():
            rank = self.rank(taxid)
            name = self.name(taxid)
            taxid = str(taxid)
            lines += '\t'.join([orf,name,rank,taxid])+'\n'
        with open(outfpath, 'w') as outfile:
            outfile.write(lines)
        if self._verbose:
            print(f'Wrote ORFs LCAs to {outfpath}')
        return outfpath

    def is_common_ancestor(self, parent, child):
        ancestor_taxid = child
        while ancestor_taxid != 1:
            if parent == ancestor_taxid:
                return True
            ancestor_taxid = self.parent(ancestor_taxid)
        return False

    def majority_vote(self, infpath, outfpath):
        """ Assigns contigs taxa using modified majority voting script
        Inputs:
            infpath - </path/to/orfs.lca>
            outfpath - </path/to/contigs.taxa>
        ReturnType: dict
        Returns:
            {contig:taxid,contig:taxid}
        """
        if not os.path.exists(infpath):
            self.write_lcas(outfpath=infpath)
        majority_vote_script = os.path.join(self.pipeline,'add_contig_taxonomy.py')
        cmd = [
            'python',
            majority_vote_script,
            self._config.get('common','ncbi_dir'),
            infpath,
            outfpath,
        ]
        if self._verbose:
            print(f'RunningMajorityVote: {" ".join(cmd)}')
        with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
            proc = subprocess.run(cmd, stdout=stdout, stderr=stderr)
        if proc.returncode:
            print(f'MajorityVoteFailed:\nArgs:{proc.args}\nReturnCode:{proc.returncode}')
            return None
        contigs = {}
        with open(outfpath) as fh:
            __ = fh.readline()
            for line in fh:
                contig,taxid = line.strip().split('\t')
                taxid = int(taxid)
                contigs.update({contig:taxid})
        self._taxa_assigned = True
        return contigs

def main(args):
    # 1. Make sure we have all of the dependencies satisfied that we need.
    dependencies.check_executables(args.verbose)
    config = dependencies.load_databases(args.config, args.verbose)
    # Likely we can use inheritance here instead of building these across
    # different classes... Especially with force,verbose,cpus,etc.
    # Should define a *base* autometa class (probably the dataset class)
    from dataset import Dataset

    dataset = Dataset(
        args.assembly,
        dbconfig=config,
        force=args.force,
        verbose=args.verbose
    )
    # Filter by length_cutoff
    dataset.filter_length(args.length_cutoff)
    # Call ORFs
    dataset.call_orfs()
    # Instatiate TaxonUtils Class
    taxutils = TaxonUtils(
        dbconfig=config,
        outdir=args.outdir,
        verbose=args.verbose,
        cpus=args.cpus,
        usepickle=args.usepickle,
    )
    # Assign Taxa to contigs
    taxutils.assign_taxa(fasta=dataset.orfs_fpath, method=args.method)
    # Write assignments to table by provided rank
    # taxutils.write_taxa(rank=args.rank)

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to retrieve taxonomy for Autometa pipeline')
    parser.add_argument('--assembly')
    parser.add_argument('--length-cutoff', default=3000, type=int)
    parser.add_argument(
        '--rank',
        default='superkingdom',
        choices=['superkingdom','phylum','class','order','family','genus','species'],
    )
    parser.add_argument('--method', default='majority_vote', choices=['majority_vote'])
    parser.add_argument('--outdir', default=os.curdir)
    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--force', action='store_true', default=False)
    parser.add_argument('--usepickle', action='store_true', default=False)
    parser.add_argument('--cpus',default=1, type=int)
    parser.add_argument(
        '--config',
        help='user defined databases config file',
        default=os.path.join(os.path.dirname(__file__), 'config', 'default.cfg'))
    args = parser.parse_args()
    main(args)
