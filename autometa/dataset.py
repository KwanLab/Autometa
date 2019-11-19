#!/usr/bin/env python3
"""
Autometa Dataset class for annotating input metagenomes
"""


import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils import GC

from taxonomy import TaxonUtils


_METHODS = {
    0:'Coverage',
    1:'ORFs',
    2:'Blastp',
    3:'LCA',
    4:'Markers',
}
NUM_METHODS = len(_METHODS)
FILTERED_EXT = 'filtered'
COV_EXT = 'cov'
ORFS_EXT = 'orfs.faa'
PRODIGAL_EXT = 'orfs.stdout'
HMMSCAN_EXT = 'hmmscan'
MARKERS_EXT = 'markers'
BLASTP_EXT = 'blastp'
LCA_EXT = 'lca'
TAXA_EXT = 'taxa'
MASTER_EXT = 'master'

FILE_EXTENSIONS = [COV_EXT, ORFS_EXT, BLASTP_EXT, LCA_EXT, MARKERS_EXT]
ORF_CALLER = 'prodigal'

class Dataset(object):
    """docstring for Autometa Dataset."""

    def __init__(self, metagenome, dbconfig, marker_set='single_copy', kingdom='bacteria', outdir=None, cpus=1, force=False, verbose=False, usepickle=False, spades=True):
        self._metagenome = os.path.realpath(metagenome)
        self._raw_dataset = os.path.realpath(metagenome)
        self._spades = spades
        self._verbose = verbose
        self._force = force
        self._usepickle = usepickle
        self._outdir = os.path.realpath(outdir) if outdir else os.path.dirname(self._metagenome)
        self._is_filtered = False
        self._seqs = SeqIO.to_dict(SeqIO.parse(self._metagenome,'fasta'))
        self._all_seqs = SeqIO.to_dict(SeqIO.parse(self._raw_dataset,'fasta'))
        self._kingdom = kingdom.lower()
        # Annotation Config & Properties
        self.cpus = str(cpus)
        self._coverage = {}
        self._all_orfs = {}
        self._lca = {}
        # Handling Marker Sets Configuration
        self._config = dbconfig
        self._marker_set = marker_set
        if not os.path.exists(self._outdir):
            if self._verbose:
                print(f'Created output directory {self._outdir}')
            os.makedirs(self._outdir)

    # ===================================
    # Definitions of Dataset Class Properties
    # ===================================
    @property
    def metagenome(self):
        return self._metagenome

    @metagenome.setter
    def metagenome(self, value):
        self._metagenome = value

    @property
    def size(self):
        return sum(len(self._seqs[rec].seq) for rec in self._seqs)

    @property
    def gc(self):
        return np.mean([GC(self._seqs[rec].seq) for rec in self._seqs])

    @property
    def coverage(self):
        return np.mean([float(seqid.split('_cov_')[-1]) for seqid in self._seqs])

    @property
    def all_seqs(self):
        return self._all_seqs

    @property
    def seqs(self):
        return self._seqs

    @seqs.setter
    def seqs(self, value):
        self._seqs = value

    @property
    def nseqs(self):
        return len(self._seqs)

    @property
    def norfs(self):
        return len(self._orfs)

    @property
    def all_orfs(self):
        return self._all_orfs

    @property
    def orfs(self):
        return self._orfs

    @orfs.setter
    def orfs(self):
        if os.path.exists(self.orfs_fpath):
            self._orfs = SeqIO.to_dict(SeqIO.parse(self.orfs_fpath, 'fasta'))
        else:
            self._orfs = {}

    @property
    def is_filtered(self):
        return self._is_filtered

    @is_filtered.setter
    def is_filtered(self, value):
        self._is_filtered = bool(value)

    def kingdom():
        doc = "The kingdom property."
        def fget(self):
            return self._kingdom
        def fset(self, value):
            self._kingdom = value
        def fdel(self):
            del self._kingdom
        return locals()
    kingdom = property(**kingdom())

    @property
    def pipeline(self):
        return os.path.join(os.path.dirname(os.path.dirname(__file__)), 'pipeline')

    # ===================================
    # Define Database Paths from Config file
    # ===================================

    @property
    def marker_set(self):
        return self._marker_set

    @marker_set.setter
    def marker_set(self, value):
        self._marker_set = value

    @property
    def hmms(self):
        hmm_config = '_'.join([self._kingdom, self.marker_set])
        return self._config.get('markers', hmm_config)

    @property
    def cutoffs(self):
        cutoffs_config = '_'.join([self._kingdom, self.marker_set, 'cutoffs'])
        return self._config.get('markers',cutoffs_config)

    @property
    def ncbi(self):
        return self._config.get('common','ncbi_dir')

    # ===================================
    # Define Output File Paths
    # ===================================

    def orfs_fpath():
        doc = "The orfs_fpath property."
        def fget(self):
            fname = '.'.join([os.path.basename(self._metagenome),ORFS_EXT])
            return os.path.join(self._outdir, fname)
        def fset(self, value):
            self._orfs_fpath = value
        def fdel(self):
            del self._orfs_fpath
        return locals()
    orfs_fpath = property(**orfs_fpath())

    def prodigal_stdout():
        doc = "The prodigal_stdout property."
        def fget(self):
            fname = '.'.join([os.path.basename(self._metagenome),PRODIGAL_EXT])
            return os.path.join(self._outdir, fname)
        def fset(self, value):
            self._prodigal_stdout = value
        def fdel(self):
            del self._prodigal_stdout
        return locals()
    prodigal_stdout = property(**prodigal_stdout())

    @property
    def hmmscan_fpath(self):
        fname = '.'.join([os.path.basename(self._metagenome),self.kingdom,HMMSCAN_EXT])
        return os.path.join(self._outdir, fname)

    @property
    def markers_fpath(self):
        fname = '.'.join([os.path.basename(self._metagenome),self.kingdom, MARKERS_EXT])
        return os.path.join(self._outdir, fname)

    @property
    def taxa_fpath(self):
        fname = '.'.join([os.path.basename(self._metagenome),TAXA_EXT])
        return os.path.join(self._outdir, fname)

    @property
    def lca_fpath(self):
        fname = '.'.join([os.path.basename(self._metagenome),LCA_EXT])
        return os.path.join(self._outdir, fname)

    @property
    def master_fpath(self):
        fname = '.'.join([os.path.basename(self._metagenome),self.kingdom,MASTER_EXT])
        return os.path.join(self._outdir, fname)


    # ===================================
    # Accessing Tables Definitions
    # ===================================

    @property
    def markers(self):
        if os.path.exists(self.markers_fpath):
            return pd.read_csv(self.markers_fpath, sep='\t', index_col='contig')
        else:
            raise FileNotFoundError(self.markers_fpath)


    @property
    def taxa(self):
        if os.path.exists(self.taxa_fpath):
            return pd.read_csv(self.taxa_fpath, sep='\t')
        else:
            raise FileNotFoundError(self.taxa_fpath)

    @property
    def lca(self):
        if os.path.exists(self.lca_fpath):
            return pd.read_csv(self.lca_fpath, sep='\t')
        else:
            raise FileNotFoundError(self.lca_fpath)

    @property
    def master(self):
        if os.path.exists(self.master_fpath):
            return pd.read_csv(self.master_fpath, sep='\t', index_col='contig')
        else:
            raise FileNotFoundError(self.master_fpath)

    # ===================================
    # Definitions of Dataset Class Methods
    # ===================================

    def filter_length(self, length_cutoff=3000):
        try:
            int(length_cutoff)
        except ValueError as err:
            print(f'{length_cutoff} can not be made into integer!')
            return
        seqs_dict = {}
        for rec in self._seqs:
            if len(self._seqs[rec].seq) >= length_cutoff:
                seqs_dict.update({self._seqs[rec].id:self._seqs[rec]})
        self._seqs = seqs_dict
        self._is_filtered = True
        fname = '.'.join([os.path.splitext(os.path.basename(self._metagenome))[0],FILTERED_EXT])
        self._metagenome = os.path.join(self._outdir, fname)

    def write_seqs(self, outfpath=''):
        if not outfpath:
            outfpath = self.metagenome
        n_written = SeqIO.write([self._seqs[rec] for rec in self._seqs], outfpath, 'fasta')
        print(f'{n_written} seqs written to {outfpath}')

    def fasta(self, records=set()):
        if not records:
            SeqIO.write([self._seqs[rec] for rec in self._seqs], sys.stdout, 'fasta')
        else:
            records = [self._seqs[rec] for rec in records if rec in self._seqs]
            SeqIO.write(records, sys.stdout, 'fasta')

    def call_orfs(self):
        # OPTIMIZE: This should use parallel or multiprocessing
        # OPTIMIZE: Piped for parallel/multiprocessing rather than writing first to file
        if not os.path.exists(self._metagenome):
            self.write_seqs()
        if os.path.exists(self.orfs_fpath) and not self._force:
            print(f'FileAlreadyExists: {self.orfs_fpath}. To overwrite use --force')
            self._orfs = SeqIO.to_dict(SeqIO.parse(self.orfs_fpath,'fasta'))
            if not self._all_orfs:
                self._all_orfs = self._orfs
            return
        cmd = ['prodigal', '-i', self._metagenome, '-a', self.orfs_fpath, '-p', 'meta', '-m', '-q']
        if self._verbose:
            print(f'RunningProdigal: {" ".join(cmd)}')
        with open(self.prodigal_stdout, 'w') as stdout, open(os.devnull, 'w') as stderr:
            proc = subprocess.run(cmd, stdout=stdout, stderr=stderr)
        if proc.returncode:
            print(f'ProdigalFailed:\nArgs:{proc.args}\nReturnCode:{proc.returncode}')
            self._orfs = {}
        else:
            self._orfs = SeqIO.to_dict(SeqIO.parse(self.orfs_fpath,'fasta'))
        if not self._all_orfs:
            self._all_orfs = self._orfs


    def get_lengths(self):
        return {seqid:len(rec.seq) for seqid,rec in self._seqs.items()}

    def get_coverages(self):
        if self._spades:
            # SPADES header: NODE_NUM_length_\d+_cov_\d+.\d+
            return {seqid:float(seqid.split('_cov_')[-1]) for seqid in self._seqs}
        else:
            raise NotImplementedError

    def get_markers(self):
        """
        Annotates markers using hmmscan on dataset ORFs and provided hmm database.
        """
        # OPTIMIZE: Similar here: Can use parallel or multiprocessing if we want to extend this
        # to grid computing (workqueue?)
        if os.path.exists(self.markers_fpath) and not self._force:
            if self._verbose:
                print(f'FileAlreadyExists: {self.markers_fpath}')
                print('force with --force if you wish to overwrite!')
            return self.markers_fpath
        cmd = ['hmmscan', '--cpu', self.cpus, '--tblout', self.hmmscan_fpath, self.hmms, self.orfs_fpath]
        if self._verbose:
            print(f'RunningHmmscan: {" ".join(cmd)}')
        with open(os.devnull, 'w') as STDOUT, open(os.devnull, 'w') as STDERR:
            proc = subprocess.run(cmd, stdout=STDOUT, stderr=STDERR)
        if proc.returncode:
            print(f'HmmscanFailed:\nArgs:{proc.args}\nReturnCode:{proc.returncode}')
            print(f'Make sure your hmm profiles are pressed!')
            print(f'hmmpress -f {self.hmms}')
            return None
        # Writes table
        # This can be optimized here but for now we will use the previous script...
        # We move from dirname=(autometa)/(dataset.py)=__file__ to ../pipeline
        markers_script = os.path.join(self.pipeline,'make_marker_table.py')
        cmd = [
            'python',
            markers_script,
            '--assembly',
            self.metagenome,
            '--hmm',
            self.hmms,
            '--cutoffs',
            self.cutoffs,
            '--processors',
            self.cpus,
            '--out',
            self.markers_fpath,
        ]
        if self._verbose:
            print(f'RunningMakeMarkerTable: {" ".join(cmd)}')
        with open(os.devnull, 'w') as STDOUT, open(os.devnull, 'w') as STDERR:
            proc = subprocess.run(cmd, stdout=STDOUT, stderr=STDERR)
        if proc.returncode:
            print(f'MakeMarkerTableFailed:\nArgs:{proc.args}\nReturnCode:{proc.returncode}')
            return None

    def filter_taxonomy(self, rank_name='bacteria', by_rank='superkingdom', method='majority_vote'):
        """ Groups contigs by provided rank using assigned taxa from method
        Inputs:
            rank_name - name within rank grouping to retrieve
            by_rank - rank to group contigs
            method - method used to assign taxonomy
        ReturnType: dict
        Returns: {rank:{contig, contig,contig,...}, rank:{...}, ...}

        Subsets ORFs and contigs to rank_name

        class variables changed:
        self._metagenome - os.path.join(self._outdir, rank_name)
        self.orfs_fpath - '.'.join([os.path.basename(self._metagenome),ORFS_EXT])
        self._seqs - subset of seqs in rank_name
        self._orfs - subset of orfs in rank_name
        """
        if self._verbose:
            print(f'Filtering Contigs by {by_rank} and taking {rank_name}')
        taxutils = TaxonUtils(
            dbconfig=self._config,
            outdir=self._outdir,
            verbose=self._verbose,
            cpus=self.cpus,
            usepickle=self._usepickle,
        )
        # taxutils.assign_taxa(fasta=self.orfs_fpath, method=method)
        contigs = taxutils.assign_taxa(
            fasta=self.orfs_fpath,
            method=method,
        )
        # Check specified rank converted from contig taxid
        rank_name = rank_name.lower()
        contig_records = []
        for contig,taxid in contigs.items():
            name = taxutils.name(taxid=taxid,rank=by_rank)
            if rank_name == name:
                contig_records.append(self._all_seqs[contig])
        if not contig_records:
            print(f'RankNotFoundError: Unable to locate {rank_name} in {by_rank}')
            return
        orf_records = []
        if not self._all_orfs:
            self.call_orfs()
        contig_ids = {record.id for record in contig_records}
        for orf,seqrecord in self._all_orfs.items():
            if orf.rsplit('_',1)[0] in contig_ids:
                orf_records.append(seqrecord)
        self._metagenome = os.path.join(self._outdir, rank_name)
        SeqIO.write(contig_records, self.metagenome, 'fasta')
        SeqIO.write(orf_records, self.orfs_fpath, 'fasta')
        self._seqs = SeqIO.to_dict(SeqIO.parse(self.metagenome,'fasta'))
        self._orfs = SeqIO.to_dict(SeqIO.parse(self.orfs_fpath,'fasta'))
        self._assigned_taxa = contigs

    def get_gc_content(self):
        return {seqid:GC(record.seq) for seqid,record in self._seqs.items()}

    def to_table(self, sep='\t'):
        if os.path.exists(self.master_fpath) and not self._force:
            print(f'FileAlreadyExists: {self.master_fpath}')
            print('To overwrite use --force')
            return self.master_fpath
        gc_content = pd.Series(self.get_gc_content())
        covs = pd.Series(self.get_coverages())
        lengths = pd.Series(self.get_lengths())
        taxids = pd.Series(self._assigned_taxa)
        df = pd.DataFrame({
            'gc':gc_content,
            'cov':covs,
            'taxid':taxids,
            'length':lengths,
        })
        df = pd.merge(
            self.markers,
            df,
            left_index=True,
            right_index=True,
            how='outer'
        )
        df.reset_index(inplace=True)
        df.to_csv(self.master_fpath, sep='\t', index=False, header=True)
        return self.master_fpath

    def get_bins(self, method='recursive_dbscan'):
        method = method.lower()
        if method == 'recursive_dbscan':
            self.recursive_dbscan()
        else:
            raise NotImplementedError

    def recursive_dbscan(self):
        script = os.path.join(self.pipeline, 'recursive_dbscan.py')
        cmd = [
            'python',
            script,
            '--input-table',
            self.master_fpath,
            '--assembly-fasta',
            self.metagenome,
            '--output-dir',
            self._outdir,
        ]
        if self._verbose:
            print(f'RunningRecursiveDBSCAN: {" ".join(cmd)}')
        with open(os.devnull, 'w') as stdout, open(os.devnull, 'w') as stderr:
            proc = subprocess.run(cmd, stdout=stdout, stderr=stderr)
        if proc.returncode:
            print(f'RecursiveDBSCANFailed:\n'
                f'Args:{proc.args}'
                f'\nReturnCode:{proc.returncode}')

if __name__ == '__main__':
    print('Dataset class for autometa pipeline... Using main.py is recommended.')
    sys.exit(0)
