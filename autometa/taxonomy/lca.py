#!/usr/bin/env python
"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

Determines the Lowest Common Ancestor given a tab-delimited BLAST table,fasta file, or iterable of SeqRecords

Assumes BLAST outfmt=6
"""


import functools
import logging
import os

import numpy as np

from tqdm import tqdm

from autometa.taxonomy.ncbi import NCBI
from autometa.common.utilities import make_pickle, unpickle, file_length
from autometa.common.external.diamond import DiamondResult
from autometa.common.external import diamond
from autometa.common.external import prodigal


logger = logging.getLogger(__name__)

class LCA(NCBI):
    """docstring for LCA."""

    def __init__(self, dbdir, outdir, usepickle=True, verbose=False, cpus=0):
        super().__init__(dbdir, verbose=verbose)
        self.outdir = outdir
        self.usepickle = usepickle
        self.verbose = verbose
        self.disable = False if verbose else True
        self.cpus = cpus
        self.tour_fp = os.path.join(self.outdir, 'tour.pkl.gz')
        self.tour = None
        self.level_fp = os.path.join(self.outdir, 'level.pkl.gz')
        self.level = None
        self.occurrence_fp = os.path.join(self.outdir, 'occurrence.pkl.gz')
        self.occurrence = None
        self.sparse_fp = os.path.join(self.outdir, 'sparse.pkl.gz')
        self.sparse = None
        self.lca_prepared = False

    def prepare_tree(self):
        """Performs Eulerian tour of nodes.dmp taxids and constructs three data
        structures:

        1. tour
        2. level
        3. occurrence

        Returns
        -------
        NoneType
            sets internals to be used for LCA lookup
        """
        if self.usepickle:
            repickle = False
            # Creating a repickle toggle b/c if *any* of the fpaths do not exist need to redo
            for fp in [self.tour_fp, self.level_fp, self.occurrence_fp]:
                if not os.path.exists(fp):
                    repickle = True
            if not repickle:
                self.tour = unpickle(fpath=self.tour_fp)
                self.level = unpickle(fpath=self.level_fp)
                self.occurrence = unpickle(fpath=self.occurrence_fp)
                return
        if self.verbose:
            logger.debug('Preparing tree, level, occurrence for LCA/RMQ')
        taxids = {}
        parents = {}
        children = {}
        # `self.nodes_fpath` inherited from NCBI
        with open(self.nodes_fpath) as fh:
            _ = fh.readline() #root
            for line in fh:
                child, parent = line.split('\t|\t')[:2]
                parent, child = [int(taxid) for taxid in [parent,child]]
                taxids.update({child:1})
                parents.update({child:parent})
                if parent in children:
                    children[parent].add(child)
                else:
                    children.update({parent:set([child])})
        tour = [(0,1)]
        direction = 1
        dist = 0
        level = [dist]
        while taxids:
            if direction > 0:
                parent = tour[-1][1]
                if parent not in children:
                    direction *= -1
                    continue
                child = children[parent].pop()
                new_node = (parent, child)
                tour.append(new_node)
                dist += 1
                level.append(dist)
                taxids.pop(child, None)
                # If the set is now empty, we need to delete the parent key in children
                if not children[parent]:
                    children.pop(parent, None)
                # Do nothing else
            elif direction < 0:
                child = tour[-1][1]
                if child not in parents:
                    direction *= -1
                    continue
                parent = parents[child]
                new_node = (child, parent)
                tour.append(new_node)
                dist -= 1
                level.append(dist)
                # Delete the child from the parents dictionary
                parents.pop(child, None)
                # If parent still has children, reverse direction
                if parent in children:
                    direction *= -1
        occurrence= {}
        for i,node in enumerate(tour):
            child = node[1]
            if child not in occurrence:
                occurrence[child]=i
        self.tour = tour
        self.level = level
        self.occurrence = occurrence
        if self.usepickle:
            make_pickle(obj=self.tour, outfpath=self.tour_fp)
            make_pickle(obj=self.level, outfpath=self.level_fp)
            make_pickle(obj=self.occurrence, outfpath=self.occurrence_fp)
        return

    def preprocess_minimums(self):
        """Preprocesses all possible LCAs and constructs into sparse table to be
         used for LCA/Range Minimum Query. Uses level array associated with its respective eulerian tour from tree construction.

        Sparse table size:
            n = number of elements in level list
            rows range = (0 to n)
            columns range = (0 to logn)

        Returns
        -------
        NoneType
            sets internals to be used for LCA lookup

        Raises
        -------
        ExceptionName
            Why the exception is raised.
        """
        if self.usepickle and os.path.exists(self.sparse_fp):
            self.sparse = unpickle(fpath=self.sparse_fp)
            return
        if self.verbose:
            logger.debug('Constructing Sparse Table')
        nrows = len(self.level)
        ncols = int(np.floor(np.log2(nrows))+1)
        sparse_matrix = np.empty((nrows,ncols))
        sparse_matrix[:,0] = [i for i in range(nrows)]
        for col in tqdm(range(1, ncols), disable=self.disable, desc='Precomputing LCAs', leave=False):
            for row in range(0, nrows):
                if 2**col > nrows:
                    continue
                if row+(2**col)-1 >= nrows:
                    sparse_matrix[row, col] = False
                    continue
                lower_index = sparse_matrix[row, (col-1)]
                upper_index = sparse_matrix[(row + 2**(col-1)), (col-1)]
                lower_index, upper_index = map(int, [lower_index, upper_index])
                lower_min = self.level[lower_index]
                upper_min = self.level[upper_index]
                if lower_min < upper_min:
                    sparse_matrix[row, col] = lower_index
                else:
                    sparse_matrix[row, col] = upper_index
        if self.usepickle:
            make_pickle(obj=sparse_matrix, outfpath=self.sparse_fp)
        self.sparse = sparse_matrix
        return

    def prepare_lca(self):
        if self.verbose:
            logger.debug('Preparing data structures for LCA')
        self.prepare_tree()
        self.preprocess_minimums()
        self.lca_prepared = True
        # tour, level, occurrence, sparse_matrix all ready
        return

    def lca(self, node1, node2):
        """Performs Range Minimum Query between 2 taxids.

        Parameters
        ----------
        node1 : int
            taxid
        node2 : int
            taxid

        Returns
        -------
        int
            LCA taxid

        Raises
        -------
        ValueError
            Provided taxid is not in the nodes.dmp tree
        """
        if not self.lca_prepared:
            self.prepare_lca()
        if node1 is None and node2 is None:
            return 1
        if node1 is None:
            return node2
        if node2 is None:
            return node1
        if node1 == node2:
            return node1
        if node1 not in self.occurrence:
            raise ValueError(f'{node1} not in tree')
        if node2 not in self.occurrence:
            raise ValueError(f'{node2} not in tree')
        if self.occurrence[node1] < self.occurrence[node2]:
            low = self.occurrence[node1]
            high = self.occurrence[node2]
        else:
            low = self.occurrence[node2]
            high = self.occurrence[node1]
        # equipartition range b/w both nodes.
        cutoff_range = int(np.floor(np.log2(high-low+1)))
        lower_index = self.sparse[low, cutoff_range]
        upper_index = self.sparse[(high-(2**cutoff_range)+1), cutoff_range]
        lower_index, upper_index = map(int, [lower_index, upper_index])
        lower_range = self.level[lower_index]
        upper_range = self.level[upper_index]
        if lower_range <= upper_range:
            lca_range = lower_range
        else:
            lca_range = upper_range
        return self.tour[self.level.index(lca_range, low, high)][1]

    def get_lcas(self, hits):
        """Retrieves the LCA of the taxids added in dict of DiamondResults

        Parameters
        ----------
        hits : dict
            {qseqid:DiamondResult, qseqid:DiamondResult, ...}
            For more details on DiamondResult see autometa.common.external.diamond

        Returns
        -------
        dict
            {qseqid:lca, qseqid:lca, ...}
        """
        lcas = {}
        n_qseqids = len(hits)
        desc = f'Determining {n_qseqids:,} qseqids\' lowest common ancestors'
        # not_found = ''
        for qseqid,hit in tqdm(hits.items(), disable=self.disable, total=n_qseqids, desc=desc, leave=False):
            taxids = set()
            for sseqid in hit.sseqids:
                taxid = hit.sseqids.get(sseqid,{'taxid':1}).get('taxid')
                if not taxid:
                    # not_found += '%s not found for %s\n' % (sseqid, qseqid)
                    # logger.warning('%s:%s suppressed/deprecated/removed', (qseqid,sseqid))
                    continue
                taxids.add(self.merged.get(taxid, taxid))
            lca = False
            num_taxids = len(taxids)
            while not lca:
                if num_taxids >= 2:
                    lca = functools.reduce(lambda taxid1,taxid2: self.lca(node1=taxid1,node2=taxid2), taxids)
                if num_taxids == 1:
                    lca = taxids.pop()
                # Exception handling where input for qseqid contains no taxids
                if num_taxids == 0:
                    lca = 1
            lcas.update({qseqid:lca})
        return lcas

    def blast2lca(self, fasta, outfpath, blast=None, hits_fpath=None, force=False):
        """Determine lowest common ancestor of provided protein sequences.

        1. Check for hits file
        2. Check for blast table
        3. Perform blast to retrieve hits make pickle for future use

        Parameters
        ----------
        fasta : str
            </path/to/prot_ORFS.fasta>.
        outfpath : str
            </path/to/lca/output/table>.
        blast : str [optional]
            </path/to/output/diamond/blastp.tsv>.
            NOTE: If provided, will assume blast has been performed and requires
            parsing (the default is None).
        hits_fpath : str [optional]
            </path/to/diamond/hits.blastp.pkl.gz>.
            hits object: {qseqid:DiamondResult, ...}.
            NOTE: If provided, will assume blast has already been performed and
            parsed and taxids added (the default is None).

        Returns
        -------
        str
            `outfpath`
        """
        if self.verbose:
            logger.debug(f'Running BLAST to LCA for {fasta}')
        if os.path.exists(outfpath) and not force:
            logger.warning(f'FileAlreadyExists {outfpath}')
            return outfpath
        if hits_fpath and os.path.exists(hits_fpath):
            hits = unpickle(hits_fpath)
            lcas = self.get_lcas(hits=hits)
            return self.write_lcas(lcas, outfpath)
        if blast and os.path.exists(blast):
            dmnd_outfpath = os.path.realpath(blast)
        else:
            fname, ext = os.path.splitext(os.path.basename(fasta))
            dmnd_outfname = '.'.join([fname,'dmnd.blastp'])
            dmnd_outfpath = os.path.join(self.outdir, dmnd_outfname)
            dmnd_outfpath = diamond.blast(
                fasta=fasta,
                database=self.nr_fpath,
                outfpath=dmnd_outfpath,
                verbose=self.verbose,
                cpus=self.cpus)
        dmnd_fname = os.path.basename(dmnd_outfpath)
        hits_fname = '.'.join([dmnd_fname, 'pkl.gz'])
        hits_fpath = os.path.join(self.outdir, hits_fname)
        hits = diamond.parse(results=dmnd_outfpath, verbose=self.verbose)
        hits = diamond.add_taxids(
            hits=hits,
            database=self.accession2taxid_fpath,
            verbose=self.verbose)
        pickled_fpath = make_pickle(obj=hits, outfpath=hits_fpath)
        lcas = self.get_lcas(hits=hits)
        return self.write_lcas(lcas, outfpath)

    def write_lcas(self, lcas, outfpath):
        """Write `lcas` to tab-delimited file: `outfpath`.
        ordered columns are: qseqid, name, rank, lca

        Parameters
        ----------
        lcas : dict
            {qseqid:lca_taxid, qseqid:lca_taxid, ...}
        outfpath : str
            </path/to/output/file.tsv>

        Returns
        -------
        str
            `outfpath`
        """
        outlines = 'qseqid\tname\trank\tlca\n'
        for qseqid,taxid in lcas.items():
            outlines += '\t'.join(map(str,[
                qseqid,
                self.name(taxid),
                self.rank(taxid),
                taxid]))+'\n'
        with open(outfpath, 'w') as outfile:
            outfile.write(outlines)
        return outfpath

    def parse(self, lca_fpath, orfs_fpath):
        """Retrieve and construct contig dictionary from provided `lca_fpath`.

        Parameters
        ----------
        lca_fpath : str
            </path/to/lcas.tsv>
            tab-delimited ordered columns: qseqid, name, rank, lca_taxid

        orfs_fpath : str
            </path/to/prodigal/called/orfs.fasta>
            Note: These ORFs should correspond to the ORFs provided in the BLAST table.

        Returns
        -------
        dict
            {contig:{rank:{taxid:counts, ...}, rank:{...}, ...}, ...}

        Raises
        -------
        FileNotFoundError
            `lca_fpath` does not exist
        FileNotFoundError
            `orfs_fpath` does not exist
        """
        logger.debug(f'Parsing LCA table: {lca_fpath}')
        if not os.path.exists(lca_fpath):
            raise FileNotFoundError(lca_fpath)
            # logger.exception(FileNotFoundError)
        if orfs_fpath and not os.path.exists(orfs_fpath):
            raise FileNotFoundError(orfs_fpath)

        translations = prodigal.contigs_from_headers(orfs_fpath)

        fname = os.path.basename(lca_fpath)
        n_lines = file_length(lca_fpath) if self.verbose else None
        disable = False if self.verbose else True
        lca_hits = {}
        with open(lca_fpath) as fh:
            header = fh.readline()
            for line in tqdm(fh, total=n_lines, disable=disable, desc=f'Parsing {fname}', leave=False):
                orf_id, name, rank, taxid = line.strip().split('\t')
                taxid = int(taxid)
                contig = translations.get(orf_id)
                if taxid != 1:
                    while rank not in set(NCBI.CANONICAL_RANKS):
                        taxid = self.parent(taxid)
                        rank = self.rank(taxid)
                if contig not in lca_hits:
                    lca_hits.update({contig:{rank:{taxid:1}}})
                elif rank not in lca_hits[contig]:
                    lca_hits[contig].update({rank:{taxid:1}})
                elif taxid not in lca_hits[contig][rank]:
                    lca_hits[contig][rank].update({taxid:1})
                else:
                    lca_hits[contig][rank][taxid] += 1
        return lca_hits

    @staticmethod
    def aggregate_lcas(self, orfs_lcas=None):
        """Aggregates ORFs' LCAS for each respective contig.

        Parameters
        ----------
        orfs_lcas : type
            Description of parameter `orfs_lcas` (the default is None).

        Returns
        -------
        dict
            {contig:{taxid:count, taxid2:count,...}, ...}

        Raises
        -------
        OrfsLcasNotFoundError
            Why the exception is raised.
        """
        logger.debug('Aggregating Contig LCAs')
        if not orfs_lcas:
            logger.exception('ORFs LCAs not found')
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

def main():
    import argparse
    basedir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    dbdir = os.path.join(basedir,'databases','ncbi')
    parser = argparse.ArgumentParser('Script to determine lowest common ancestor')
    parser.add_argument('orfs', help='<path/to/orfs.fasta>')
    parser.add_argument('blast', help='<path/to/blast.tsv>')
    parser.add_argument('outdir', help='<path/to/output/dir>')
    parser.add_argument('outfname', help='<lca filename>')
    parser.add_argument('--blast-hits', help='<path/to/blast.pkl.gz> (with taxids already added)')
    parser.add_argument('--dbdir', help='<path/to/ncbi/dir>', default=dbdir)
    parser.add_argument('--nopickle', help='do not pickle objects to disk',
        action='store_false', default=True)
    parser.add_argument('--verbose', help="add verbosity", action='store_true', default=False)
    parser.add_argument('--force', help="force overwrite if file already exists", action='store_true', default=False)
    args = parser.parse_args()
    lca = LCA(
        dbdir=args.dbdir,
        outdir=args.outdir,
        usepickle=args.nopickle,
        verbose=args.verbose,
    )
    outfpath = os.path.join(args.outdir, args.outfname)
    lca.blast2lca(
        fasta=args.orfs,
        outfpath=outfpath,
        blast=args.blast,
        hits_fpath=args.blast_hits,
        force=args.force)

if __name__ == '__main__':
    main()
