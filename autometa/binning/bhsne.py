#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Run the python2 implementation of tSNE for embedding comparisons.

SETUP ENV:
conda create -n py2 python=2.7
conda activate py2
conda install -c anaconda scikit-learn atlas cython -y
conda update -c anaconda scikit-learn -y
conda install -c maxibor tsne -y
conda install pandas -y
"""


import logging

import pandas as pd

from sklearn.decomposition import PCA
from tsne import bh_sne


logger = logging.getLogger(__name__)


def embed(kmers_fpath):
    """Embed k-mers using python2.7 implementation of tSNE

    Parameters
    ----------
    kmer_fpath : str
        </path/to/kmers.normalized.tsv>

    Returns
    -------
    pd.DataFrame
        embedded kmers

    Raises
    -------
    FileNotFoundError
        `kmers_fpath` does not exist

    """
    df = pd.read_csv(kmers_fpath, sep='\t', index_col='contig')
    df.dropna(how='all',inplace=True)
    X = df.to_numpy()
    X = PCA(n_components=50).fit_transform(X)
    X = bh_sne(X, d=2)
    return pd.DataFrame(X, columns=['x','y'], index=df.index)

def main(args):
    df = embed(args.kmers)
    logger.debug('{} embedded. : df.shape: {}'.format(args.kmers, df.shape))
    df.to_csv(args.embedded, sep='\t', index=True, header=True)
    logger.debug('embedded written {}'.format(args.embedded))

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('py2.7 tsne manifold learning implementation for k-mer frequency embedding.')
    parser.add_argument('kmers',help='</path/to/kmers.normalized.tsv>')
    parser.add_argument('embedded',help='</path/to/kmers.embedded.tsv>')
    args = parser.parse_args()
    main(args)
