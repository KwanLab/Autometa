#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script containing Metagenome class for general handling of assembly
"""


import os

from Bio import SeqIO

from .mag import Mag

class Metagenome:
    """docstring for Metagenome."""

    def __init__(self, assembly):
        self.assembly = assembly
        self.dirname = os.path.dirname(os.path.realpath(self.assembly))
        self.sequences = [seq for seq in SeqIO.parse(self.assembly, 'fasta')]
        self.orfs_called = False
        self.kmers = None
        self.taxonomy = None

    def __str__(self):
        return self.assembly

    @property
    def nucls(self):
        return self.get_orfs(orf_type='nucl')

    @property
    def prots(self):
        return self.get_orfs(orf_type='prot')

    def length_filter(self, cutoff=3000):
        """Filters sequences by length with provided cutoff
        Input:
            cutoff - int|float g.t. 0
        ReturnType: generator
        Returns: generator object of sequnces with cutoff filter applied
        """
        if not type(cutoff) in [int, float]:
            raise TypeError(f'{cutoff} must be a float or int')
        if not cutoff > 0:
            raise ValueError(f'cutoff must be a positive real number')
        return Mag(
            parent_metagenome=self.assembly,
            contig_names=[seq.id for seq in self.sequences if len(seq) >= cutoff])

    def call_orfs(self, force=False, verbose=True):
        """Calls ORFs using external executable: Prodigal
        Inputs:
            force - force overwrite of existing ORFs files (Boolean)
            verbose - add verbosity (Boolean)
        ReturnType: 3-tuple
        Returns:
            on success - (True, nucls_filepath, prots_filepath)
            on failure - (False, None, None)
        """
        for arg in [force, verbose]:
            if type(arg) is not bool:
                raise TypeError(f'{arg} must be a boolean. I.e. True|False')

        from .common.external.prodigal import get_orfs
        # OPTIMIZE: Should not need to call ORFs on contigs below length cutoff
        self.basename = os.path.basename(self.assembly)
        nucls_ext = 'orfs.fna'
        prots_ext = 'orfs.faa'
        split_bname = os.path.splitext(self.basename)[0]
        self.nucls_fname = '.'.join([split_bname, nucls_ext])
        self.prots_fname = '.'.join([split_bname, prots_ext])
        self.nucls_out = os.path.join(self.dirname, self.nucls_fname)
        self.prots_out = os.path.join(self.dirname, self.prots_fname)
        success, nucls_fp, prots_fp = get_orfs(
            assembly=self.assembly,
            nucls_out=self.nucls_out,
            prots_out=self.prots_out,
            force=force,
            verbose=verbose)
        if not success:
            print('Calling ORFs with Prodigal Failed!')
        else:
            self.orfs_called = True

    def get_orfs(self, orf_type='prot'):
        """Gets ORFs after being called from self.call_orfs
        Inputs:
            orf_type - format of ORFs to retrieve choices=['nucl','prot']
                either nucleotide or amino acids
        ReturnType: List
        Returns [SeqRecord, ...]
        """
        if not self.orfs_called:
            self.call_orfs()
        if orf_type not in {'prot','nucl'}:
            raise ValueError('orf_type must be prot|nucl!')
        orfs_fpath = self.prots if orf_type == 'prot' else self.nucls
        return [orf for orf in SeqIO.parse(orfs_fpath, 'fasta')]

    def get_kmers(self, ):
        raise NotImplementedError

    def get_taxonomy(self, ):
        raise NotImplementedError

def main(args):
    mg = Metagenome(args.assembly)
    mag = mg.length_filter(cutoff=args.cutoff)
    print(
        f'{mag.parent_metagenome} -> filtered by cutoff: {args.cutoff}\n'
        f'{mag.n_contigs} seqs remain after filter cutoff'
        f'Contig Names:\n{mag.contig_names}')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Metagenome class to filter sequences by length')
    parser.add_argument('assembly', help='</path/to/assembly.fasta>')
    parser.add_argument('--cutoff', help='length to filter sequences',default=3000)
    args = parser.parse_args()
    main(args)
