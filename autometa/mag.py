#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Autometa Bin Class
"""

import numpy as np

from Bio import SeqIO


class Mag:
    """docstring for Mag."""

    def __init__(self, parent_metagenome, contig_names):
        self.parent_metagenome = parent_metagenome
        self.contig_names = contig_names
        self.marker_set = None
        self.core_markers = None
        self.is_finished = False

    def __str__(self):
        return ','.join(self.contig_names)

    # def __len__(self):
    #     return np.sum(len(seq) for seq in self.get_seqs())

    @property
    def n_contigs(self):
        return len(self.contig_names)

    @property
    def contigs(self):
        return self.get_seqs()

    def get_seqs(self):
        return [seq for seq in SeqIO.parse(self.parent_metagenome, 'fasta')
            if seq.id in self.contig_names]

    def split_taxonomy(self, ):
        raise NotImplementedError

    def split_nucleotides(self, ):
        raise NotImplementedError


def main(args):
    print('Bin class used by Autometa Metagenomics Binning Pipeline')
    pass

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser('Autometa Bin Class')
    parser.add_argument('--test')
    args = parser.parse_args()
    main(args)
