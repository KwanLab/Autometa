#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File containing functions to count, retrieve, k-mers given sequences

TODO: Separate file to handle parallel,work-queue processing
"""


import logging
import os

import numpy as np
import pandas as pd
import multiprocessing as mp

from tqdm import tqdm
from Bio import SeqIO
from scipy.stats import gmean
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from umap import UMAP

from autometa.common.utilities import gunzip

logger = logging.getLogger(__name__)
# TODO: Silence UMAP DEBUG logger (numba)
# logger = logging.getLogger(UMAP.__name__).setLevel(logging.WARNING)

def revcomp(string):
    complement = {'A':'T','T':'A','C':'G','G':'C'}
    complements = []
    for i in range(len(string)):
        if string[i] in complement:
            complements.append(complement[string[i]])
        else:
            return -1
    return ''.join(reversed(complements))

def init_kmers(kmer_size=5):
    # Count K-mer frequencies
    # Holds lists of k-mer counts, keyed by contig name
    kmers = {}
    counts = 0
    uniq_kmers = dict()
    dna_letters = ['A', 'T', 'C', 'G']
    all_kmers = list(dna_letters)
    for i in range(1, kmer_size):
        new_list = []
        for current_seq in all_kmers:
            for char in dna_letters:
                new_list.append(current_seq + char)
        all_kmers = new_list
    # Now we trim k-mers and put them in the dictionary
    # Q: What is being trimmed here?
    # A: I think trim means subset by unique from the reverse complement
    for kmer in all_kmers:
        kmer_reverse = revcomp(kmer)
        if (kmer not in uniq_kmers) and (kmer_reverse not in uniq_kmers):
            uniq_kmers[kmer] = counts
            counts += 1
    return uniq_kmers

def load(kmers_fpath):
    """Load in a previously counted k-mer frequencies table.

    Parameters
    ----------
    kmers_fpath : str
        Description of parameter `kmers_fpath`.

    Returns
    -------
    pd.DataFrame
        index='contig', cols=[kmer, kmer, ...]

    Raises
    -------
    FileNotFoundError
        `kmers_fpath` does not exist or is empty

    """
    if not os.path.exists(kmers_fpath) or os.stat(kmers_fpath).st_size == 0:
        raise FileNotFoundError(kmers_fpath)
    return pd.read_csv(kmers_fpath, sep='\t', index_col='contig')

def mp_counter(assembly, ref_kmers, nproc=mp.cpu_count()):
    pool = mp.Pool(nproc)
    args = [(record,ref_kmers) for record in SeqIO.parse(assembly, 'fasta')]
    logger.debug(f'Pool counter (nproc={nproc}): counting {len(args):,} records k-mer frequencies')
    results = pool.map(kmer_counter, args)
    pool.close()
    pool.join()
    return results

def kmer_counter(args):
    record, ref_kmers = args
    for ref_kmer in ref_kmers:
        kmer_size = len(ref_kmer)
        break
    n_uniq_kmers = len(ref_kmers)
    contig_kmer_counts = [0] * n_uniq_kmers
    record_length = len(record.seq)
    max_length = record_length - kmer_size
    if max_length <= 0:
        logger.warning(f'{record.id} can not be counted! k-mer size exceeds length. {record_length}')
        return {record.id:contig_kmer_counts}
        # contig_kmer_counts.insert(0, record.id)
        # return '\t'.join([str(c) for c in contig_kmer_counts])+'\n'
    for i in range(max_length):
        kmer = record.seq[i:i+kmer_size]
        kmer_revcomp = kmer.reverse_complement()
        kmer, kmer_revcomp = map(str, [kmer,kmer_revcomp])
        if kmer not in ref_kmers and kmer_revcomp not in ref_kmers:
            continue
        if kmer in ref_kmers:
            index = ref_kmers[kmer]
        else:
            index = ref_kmers[kmer_revcomp]
        contig_kmer_counts[index] += 1
    # contig_kmer_counts.insert(0, record.id)
    return {record.id:contig_kmer_counts}

def seq_counter(assembly, ref_kmers, verbose=True):
    n_uniq_kmers = len(ref_kmers)
    for ref_kmer in ref_kmers:
        kmer_size = len(ref_kmer)
        break
    kmer_counts = {}
    seqs = SeqIO.parse(assembly, 'fasta')
    desc = f'Counting {n_uniq_kmers} unique {kmer_size}-mers'
    disable = not verbose
    for record in tqdm(seqs, desc=desc, disable=disable, leave=False):
        contig_kmer_counts = [0] * n_uniq_kmers
        max_length = len(record.seq) - kmer_size
        if max_length <= 0:
            logger.warning(f'{record.id} can not be counted! k-mer size exceeds length. {len(record.seq)}')
            kmer_counts.update({record.id:contig_kmer_counts})
            continue
        for i in range(max_length):
            kmer = record.seq[i:i+kmer_size]
            kmer_revcomp = kmer.reverse_complement()
            kmer, kmer_revcomp = map(str, [kmer,kmer_revcomp])
            if kmer not in ref_kmers and kmer_revcomp not in ref_kmers:
                continue
            if kmer in ref_kmers:
                index = ref_kmers[kmer]
            else:
                index = ref_kmers[kmer_revcomp]
            contig_kmer_counts[index] += 1
        kmer_counts.update({record.id:contig_kmer_counts})
    return kmer_counts

def count(assembly, kmer_size=5, normalized=False, verbose=True, multiprocess=True,
    nproc=mp.cpu_count()):
    """Counts k-mer frequencies for provided assembly file

    First we make a dictionary of all the possible k-mers (discounting reverse
    complements). Each k-mer's count is updated by index when encountered in the
    record.

    Parameters
    ----------
    assembly : str
        Description of parameter `assembly`.
    kmer_size : int
        length of k-mer to count `kmer_size` (the default is 5).
    normalized : bool
        Whether to return the CLR normalized dataframe (the default is True).

    Returns
    -------
    pandas.DataFrames
        index_col='contig', tab-delimited, cols=unique_kmers
        i.e. 5-mer columns=[AAAAA, AAAAT, AAAAC, AAAAG,  ..., GCCGC]

    Raises
    -------
    TypeError
        `kmer_size` must be an int

    """
    if not type(kmer_size) is int:
        raise TypeError(f'kmer_size must be an int! Given: {kmer_size}')
    ref_kmers = init_kmers(kmer_size)
    if assembly.endswith('.gz'):
        assembly = gunzip(assembly, assembly.rstrip('.gz'))
    if multiprocess:
        kmer_counts = {}
        counts = mp_counter(assembly, ref_kmers, nproc)
        for result in counts:
            kmer_counts.update(result)
    else:
        kmer_counts = seq_counter(assembly, ref_kmers, verbose=verbose)
    df = pd.DataFrame(kmer_counts, index=ref_kmers).transpose()
    df.index.name = 'contig'
    if normalized:
        return normalize(df)
    else:
        return df

def normalize(df):
    """Normalize k-mers by Centered Log Ratio transformation

    1. Drop any k-mers not contained by any contigs
    2a. Normalize the k-mer count by the total count of all k-mers for a given contig
    2b. Add 1 as 0 can not be utilized for CLR
    3. Perform CLR transformation log(norm. value / geometric mean norm. value)

    Parameters
    ----------
    df : pd.DataFrame
        K-mers Dataframe where index_col='contig' and column values are k-mer
        frequencies.

    References:
    - Aitchison, J. The Statistical Analysis of Compositional Data (1986)
    - Pawlowsky-Glahn, Egozcue, Tolosana-Delgado. Lecture Notes on Compositional Data Analysis (2011)
    - https://stats.stackexchange.com/questions/242445/why-is-isometric-log-ratio-transformation-preferred-over-the-additivealr-or-ce
    - https://stats.stackexchange.com/questions/305965/can-i-use-the-clr-centered-log-ratio-transformation-to-prepare-data-for-pca
    - http://www.sediment.uni-goettingen.de/staff/tolosana/extra/CoDa.pdf

    Returns
    -------
    pd.DataFrame
        index='contig', cols=[kmer, kmer, ...]
        Columns have been transformed by CLR normalization.
    """
    # OPTIMIZE: May be able to implement this transformation with dask?
    return df.dropna(axis='columns', how='all')\
        .transform(lambda x: (x+1) / x.sum(), axis='columns')\
        .transform(lambda x: np.log(x / gmean(x)), axis='columns')


def embed(kmers=None, embedded=None, n_components=3, do_pca=True, pca_dimensions=50, method='TSNE', perplexity=30, **kwargs):
    """Embed k-mers using provided `method`.

    Parameters
    ----------
    kmers : str or pd.DataFrame
        Description of parameter `kmers` (the default is None).
    embedded : str
        Description of parameter `embedded` (the default is None).
    n_components : int
        `n_components` to embed k-mer frequencies (the default is 3).
    do_pca : bool
        Perform PCA decomposition prior to embedding (the default is True).
    pca_dimensions : int
        Reduce k-mer frequencies dimensions to `pca_dimensions` (the default is 50).
        If None, will estimate based on
    method : str
        Description of parameter `method` (the default is 'TSNE').
    perplexity : float
        Description of parameter `perplexity` (the default is 30).
    **kwargs : dict
        Other keyword arguments to be supplied to respective `method`.

    Returns
    -------
    pd.DataFrame
        embedded dataframe with index='contig' and cols=['x','y','z']

    Raises
    -------
    ValueError
        Either `kmers` or `embedded` must be provided.
    ValueError
        Provided `method` is not an available choice.
    FileNotFoundError
        `kmers` type must be a pd.DataFrame or filepath.
    """
    if not kmers and not embedded:
        raise ValueError('kmers or embedded is required')
    df = None
    if kmers and type(kmers) is str and os.path.exists(kmers) and os.stat(kmers).st_size >0:
        df = pd.read_csv(kmers, sep='\t', index_col='contig')
    elif kmers and type(kmers) is pd.DataFrame:
        df = kmers
    if embedded and os.path.exists(embedded) and os.stat(embedded).st_size > 0:
        logger.debug(f'k-mers frequency embedding already exists {embedded}')
        return pd.read_csv(embedded, sep='\t', index_col='contig')
    if df is None or df.empty:
        kmers_desc = f'kmers:{kmers} type:{type(kmers)}'
        embed_desc = f'embedded:{embedded} type:{type(embedded)}'
        requirements = f'kmers type must be a pd.DataFrame or filepath.'
        raise FileNotFoundError(f'{kmers_desc} {embed_desc} {requirements}')

    method = method.upper()
    if method not in ['UMAP','TSNE']:
        raise ValueError(f'{method} not in embedding methods. Choices: TSNE, UMAP')
    # PCA
    n_samples, n_dims = df.shape
    # Drop any rows that all cols contain NaN. This may occur if the contig length is below the k-mer size
    df.dropna(how='all', inplace=True)
    X = df.to_numpy()
    if n_dims > pca_dimensions and do_pca:
        logger.debug(f'Performing decomposition with PCA: {n_dims} to {pca_dimensions} dims')
        X = PCA(n_components=pca_dimensions).fit_transform(X)
        # X = PCA(n_components='mle').fit_transform(X)
        n_samples, n_dims = X.shape

    logger.debug(f'{method}: {n_samples} data points and {n_dims} dimensions')
    # Adjust perplexity according to the number of data points
    n_rows = n_samples-1
    scaler = 3.0
    if n_rows < (scaler*perplexity):
        perplexity = (n_rows/scaler) - 1

    def do_TSNE(perplexity=perplexity, n_components=n_components):
        return TSNE(
            n_components=n_components,
            perplexity=perplexity,
            random_state=0).fit_transform(X)

    def do_UMAP(n_neighbors=15, n_components=n_components, metric='euclidean'):
        return UMAP(
            n_neighbors=n_neighbors,
            n_components=n_components,
            metric=metric).fit_transform(X)

    dispatcher = {'TSNE':do_TSNE, 'UMAP':do_UMAP}
    logger.debug(f'Performing embedding with {method}')
    X = dispatcher[method](**kwargs)
    if n_components == 3:
        embedded_df = pd.DataFrame(X, columns=['x','y','z'], index=df.index)
    elif n_components == 2:
        embedded_df = pd.DataFrame(X, columns=['x','y'], index=df.index)
    else:
        embedded_df = pd.DataFrame(X, index=df.index)
    if embedded:
        embedded_df.to_csv(embedded, sep='\t', index=True, header=True)
        logger.debug(f'embedded.shape {embedded_df.shape} : Written {embedded}')
    return embedded_df

def main(args):
    try:
        df = load(args.kmers)
        logger.debug(f'{args.kmers} exists... loaded: df.shape {df.shape}')
    except FileNotFoundError as err:
        df = count(
            assembly=args.fasta,
            kmer_size=args.size,
            normalized=False,
            multiprocess=args.multiprocess,
            nproc=args.nproc)
        df.to_csv(args.kmers, sep='\t', header=True, index=True)
        logger.debug(f'Wrote {len(df)} contigs {args.size}-mers frequencies to {args.kmers}.')

    if args.normalized:
        try:
            ndf = load(args.normalized)
            logger.debug(f'{args.normalized} exists... loaded: df.shape {ndf.shape}')
        except FileNotFoundError as err:
            logger.debug(f'Normalizing {df.shape} k-mers DataFrame.')
            ndf = normalize(df)
            ndf.to_csv(args.normalized, sep='\t', header=True, index=True)
            logger.debug(f'Wrote {len(df)} normalized k-mer freqs. to {args.normalized}.')

    if not args.embedded:
        import sys;sys.exit(0)

    if args.normalized:
        logger.debug(f'Embedding {args.normalized}')
        embedded_df = embed(
            kmers=args.normalized,
            embedded=args.embedded,
            method=args.method)
    else:
        logger.debug(f'Embedding {args.kmers}')
        embedded_df = embed(
            kmers=args.kmers,
            embedded=args.embedded,
            method=args.method)

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Count k-mers')
    parser.add_argument('--fasta', help='</path/to/sequences.fna>', required=True)
    parser.add_argument('--kmers', help='</path/to/kmers.tsv>', required=True)
    parser.add_argument('--size', help='k-mer size', default=5, type=int)
    parser.add_argument('--normalized', help='</path/to/kmers.normalized.tsv>')
    parser.add_argument('--embedded', help='</path/to/kmers.embedded.tsv>')
    parser.add_argument('--method', help='embedding method', choices=['TSNE','UMAP'], default='UMAP')
    parser.add_argument('--multiprocess', help='count k-mers using multiprocessing',
        action='store_true', default=False)
    parser.add_argument('--nproc', help='num. processors to use if multiprocess is selected',
        default=mp.cpu_count(), type=int)
    args = parser.parse_args()
    main(args)
