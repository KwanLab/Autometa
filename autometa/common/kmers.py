#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

Count, normalize and embed k-mers given nucleotide sequences
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
from tsne import bh_sne
from umap import UMAP

from autometa.common.utilities import gunzip
from autometa.common.exceptions import KmerFormatError
from autometa.common.exceptions import KmerEmbeddingError

# Suppress numba logger debug output
numba_logger = logging.getLogger("numba").setLevel(logging.ERROR)
logger = logging.getLogger(__name__)


def _revcomp(string):
    """Reverse complement the provided `string`.

    Parameters
    ----------
    string : str
        A k-mer string generated from `init_kmers`

    Returns
    -------
    str
        reverse complemented string.

    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement.get(char) for char in reversed(string))


def init_kmers(kmer_size=5):
    """Initialize k-mers from `kmer_size`. Respective reverse complements will
    be removed.

    Parameters
    ----------
    kmer_size : int, optional
        pattern size of k-mer to intialize dict (the default is 5).

    Returns
    -------
    dict
        {kmer:index, ...}

    """
    kmers = {}
    index = 0
    uniq_kmers = dict()
    dna_letters = ["A", "T", "C", "G"]
    all_kmers = list(dna_letters)
    for i in range(1, kmer_size):
        new_list = []
        for current_seq in all_kmers:
            for char in dna_letters:
                new_list.append(current_seq + char)
        all_kmers = new_list
    # subset uniq_kmers by removing any reverse complements
    for kmer in all_kmers:
        kmer_reverse = _revcomp(kmer)
        if (kmer not in uniq_kmers) and (kmer_reverse not in uniq_kmers):
            uniq_kmers[kmer] = index
            index += 1
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
    KmerFormatError
        `kmers_fpath` file format is invalid

    """
    if not os.path.exists(kmers_fpath) or os.stat(kmers_fpath).st_size == 0:
        raise FileNotFoundError(kmers_fpath)
    try:
        df = pd.read_csv(kmers_fpath, sep="\t", index_col="contig")
    except ValueError as err:
        raise KmerFormatError(kmers_fpath) from ValueError
    return df


def mp_counter(assembly, ref_kmers, nproc=mp.cpu_count()):
    """Multiprocessing k-mer counter used in `count`. (Should not be used directly).

    Parameters
    ----------
    assembly : str
        </path/to/assembly.fasta> (nucleotides)
    ref_kmers : dict
        {kmer:index, ...}
    nproc : int, optional
        Number of cpus to use. (the default will use all available).

    Returns
    -------
    list
        [{record:counts}, {record:counts}, ...]

    """
    pool = mp.Pool(nproc)
    args = [(record, ref_kmers) for record in SeqIO.parse(assembly, "fasta")]
    logger.debug(
        f"Pool counter (nproc={nproc}): counting {len(args):,} records k-mer frequencies"
    )
    results = pool.map(record_counter, args)
    pool.close()
    pool.join()
    return results


def record_counter(args):
    """single record counter used when multiprocessing.

    Parameters
    ----------
    args : 2-tuple
        (record, ref_kmers)
        - record : SeqIO.SeqRecord
        - ref_kmers : {kmer:index, ...}

    Returns
    -------
    dict
        {contig:[count,count,...]} count index is respective to ref_kmers.keys()

    """
    record, ref_kmers = args
    for ref_kmer in ref_kmers:
        kmer_size = len(ref_kmer)
        break
    n_uniq_kmers = len(ref_kmers)
    contig_kmer_counts = [0] * n_uniq_kmers
    record_length = len(record.seq)
    max_length = record_length - kmer_size
    if max_length <= 0:
        logger.warning(
            f"{record.id} can not be counted! k-mer size exceeds length. {record_length}"
        )
        contig_kmer_counts = [pd.NA] * n_uniq_kmers
        return {record.id: contig_kmer_counts}
    for i in range(max_length):
        kmer = record.seq[i : i + kmer_size]
        # reverse_complement() is Biopython specific method for SeqRecord object
        kmer_revcomp = kmer.reverse_complement()
        kmer, kmer_revcomp = map(str, [kmer, kmer_revcomp])
        if kmer not in ref_kmers and kmer_revcomp not in ref_kmers:
            continue
        if kmer in ref_kmers:
            index = ref_kmers[kmer]
        else:
            index = ref_kmers[kmer_revcomp]
        contig_kmer_counts[index] += 1
    contig_kmer_counts = [c if c != 0 else pd.NA for c in contig_kmer_counts]
    return {record.id: contig_kmer_counts}


def seq_counter(assembly, ref_kmers, verbose=True):
    """Sequentially count k-mer frequencies.

    Parameters
    ----------
    assembly : str
        </path/to/assembly.fasta> (nucleotides)
    ref_kmers : dict
        {kmer:index, ...}
    verbose : bool, optional
        enable progress bar `verbose` (the default is True).

    Returns
    -------
    dict
        {contig:[count,count,...]} count index is respective to ref_kmers.keys()

    """
    n_uniq_kmers = len(ref_kmers)
    for ref_kmer in ref_kmers:
        kmer_size = len(ref_kmer)
        break
    kmer_counts = {}
    seqs = SeqIO.parse(assembly, "fasta")
    desc = f"Counting {n_uniq_kmers} unique {kmer_size}-mers"
    disable = not verbose
    for record in tqdm(seqs, desc=desc, disable=disable, leave=False):
        contig_kmer_counts = [0] * n_uniq_kmers
        max_length = len(record.seq) - kmer_size
        if max_length <= 0:
            logger.warning(
                f"{record.id} can not be counted! k-mer size exceeds length. {len(record.seq)}"
            )
            contig_kmer_counts = [pd.NA] * n_uniq_kmers
            kmer_counts.update({record.id: contig_kmer_counts})
            continue
        for i in range(max_length):
            kmer = record.seq[i : i + kmer_size]
            # reverse_complement() is Biopython specific method for SeqRecord object
            kmer_revcomp = kmer.reverse_complement()
            kmer, kmer_revcomp = map(str, [kmer, kmer_revcomp])
            if kmer not in ref_kmers and kmer_revcomp not in ref_kmers:
                continue
            if kmer in ref_kmers:
                index = ref_kmers[kmer]
            else:
                index = ref_kmers[kmer_revcomp]
            contig_kmer_counts[index] += 1
        contig_kmer_counts = [c if c != 0 else pd.NA for c in contig_kmer_counts]
        kmer_counts.update({record.id: contig_kmer_counts})
    return kmer_counts


def count(
    assembly,
    kmer_size=5,
    normalized=False,
    verbose=True,
    multiprocess=True,
    nproc=mp.cpu_count(),
):
    """Counts k-mer frequencies for provided assembly file

    First we make a dictionary of all the possible k-mers (discounting reverse
    complements). Each k-mer's count is updated by index when encountered in the
    record.

    Parameters
    ----------
    assembly : str
        Description of parameter `assembly`.
    kmer_size : int, optional
        length of k-mer to count `kmer_size` (the default is 5).
    normalized : bool, optional
        Whether to return the CLR normalized dataframe (the default is True).
    verbose : bool, optional
        Enable progress bar `verbose` (the default is True).
    multiprocess : bool, optional
        Use multiple cores to count k-mer frequencies (the default is True).
    nproc : int, optional
        Number of cpus to use. (the default will use all available).

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
        raise TypeError(f"kmer_size must be an int! Given: {type(kmer_size)}")
    ref_kmers = init_kmers(kmer_size)
    if assembly.endswith(".gz"):
        assembly = gunzip(assembly, assembly.rstrip(".gz"))
    if multiprocess:
        kmer_counts = {}
        counts = mp_counter(assembly, ref_kmers, nproc)
        for result in counts:
            kmer_counts.update(result)
    else:
        kmer_counts = seq_counter(assembly, ref_kmers, verbose=verbose)
    df = pd.DataFrame(kmer_counts, index=ref_kmers).transpose()
    df.index.name = "contig"
    if normalized:
        return normalize(df)
    else:
        return df


def normalize(df):
    """Normalize k-mers by Centered Log Ratio transformation

    1a. Drop any k-mers not present for all contigs
    1b. Drop any contigs not containing any kmer counts
    1c. Fill any remaining na values with 0
    2a. Normalize the k-mer count by the total count of all k-mers for a given contig
    2b. Add 1 as 0 can not be utilized for CLR
    3. Perform CLR transformation log(norm. value / geometric mean norm. value)

    Parameters
    ----------
    df : pd.DataFrame
        K-mers Dataframe where index_col='contig' and column values are k-mer
        frequencies.

    # TODO: Place these references in readthedocs documentation and remove from def.
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
    # steps in 1: data cleaning
    df.dropna(axis="columns", how="all", inplace=True)
    df.dropna(axis="index", how="all", inplace=True)
    df.fillna(0, inplace=True)
    # steps in 2 and 3: normalization and CLR transformation
    step_2a = lambda x: (x + 1) / x.sum()
    step_2b = lambda x: np.log(x / gmean(x))
    return df.transform(step_2a, axis="columns").transform(step_2b, axis="columns")


def embed(
    kmers=None,
    embedded=None,
    n_components=2,
    do_pca=True,
    pca_dimensions=50,
    method="bhsne",
    perplexity=30.0,
    **kwargs,
):
    """Embed k-mers using provided `method`.

    Notes
    -----

        * `sklearn.manifold.TSNE <https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html#sklearn.manifold.TSNE>`_
        * `UMAP <https://umap-learn.readthedocs.io/en/latest/>`_
        * `tsne.bh_sne <https://pypi.org/project/tsne/>`_

    Parameters
    ----------
    kmers : str or pd.DataFrame
        </path/to/input/kmers.normalized.tsv>
    embedded : str, optional
        </path/to/output/kmers.embedded.tsv> If provided will write to `embedded`.
    n_components : int, optional
        `n_components` to embed k-mer frequencies (the default is 2).
    do_pca : bool, optional
        Perform PCA decomposition prior to embedding (the default is True).
    pca_dimensions : int, optional
        Reduce k-mer frequencies dimensions to `pca_dimensions` (the default is 50).
        If None, will estimate based on
    method : str, optional
        embedding method to use (the default is 'bhsne').
        choices include sksne, bhsne and umap.
    perplexity : float, optional
        hyperparameter used to tune sksne and bhsne (the default is 30.0).
    **kwargs : dict, optional
        Other keyword arguments to be supplied to respective `method`.

    Returns
    -------
    pd.DataFrame
        embedded dataframe with index='contig' and cols=['x','y','z']

    Raises
    -------
    KmerEmbeddingError
        Either `kmers` or `embedded` must be provided.
    KmerFormatError
        Provided `kmers` or `embedded` are not formatted correctly for use.
    ValueError
        Provided `method` is not an available choice.
    FileNotFoundError
        `kmers` type must be a pd.DataFrame or filepath.
    """
    if not kmers and not embedded:
        msg = f"`kmers` (given: {kmers}) or `embedded` (given: {embedded}) is required"
        raise KmerEmbeddingError(msg)
    df = None
    if (
        kmers
        and type(kmers) is str
        and os.path.exists(kmers)
        and os.stat(kmers).st_size > 0
    ):
        try:
            df = pd.read_csv(kmers, sep="\t", index_col="contig")
        except ValueError as err:
            raise KmerFormatError(embedded) from ValueError
    elif kmers and type(kmers) is pd.DataFrame:
        df = kmers
    if embedded and os.path.exists(embedded) and os.stat(embedded).st_size > 0:
        logger.debug(f"k-mers frequency embedding already exists {embedded}")
        try:
            df = pd.read_csv(embedded, sep="\t", index_col="contig")
        except ValueError as err:
            raise KmerFormatError(embedded) from ValueError
        return df

    if df is None or df.empty:
        kmers_desc = f"kmers:{kmers} type:{type(kmers)}"
        embed_desc = f"embedded:{embedded} type:{type(embedded)}"
        requirements = f"kmers type must be a pd.DataFrame or filepath."
        raise FileNotFoundError(f"{kmers_desc} {embed_desc} {requirements}")

    method = method.lower()
    choices = {"umap", "sksne", "bhsne"}
    if method not in choices:
        raise ValueError(
            f'{method} not in embedding methods. Choices: {", ".join(choices)}'
        )
    # PCA
    n_samples, n_dims = df.shape
    # Drop any rows that all cols contain NaN. This may occur if the contig length is below the k-mer size
    df.dropna(axis="index", how="all", inplace=True)
    df.fillna(0, inplace=True)
    X = df.to_numpy()
    if n_dims > pca_dimensions and do_pca:
        logger.debug(
            f"Performing decomposition with PCA: {n_dims} to {pca_dimensions} dims"
        )
        X = PCA(n_components=pca_dimensions).fit_transform(X)
        # X = PCA(n_components='mle').fit_transform(X)
        n_samples, n_dims = X.shape

    logger.debug(f"{method}: {n_samples} data points and {n_dims} dimensions")

    n_rows = n_samples - 1
    scaler = 3.0
    if n_rows < (scaler * perplexity):
        perplexity = (n_rows / scaler) - 1

    def do_sksne(perplexity=perplexity, n_components=n_components, seed=42):
        # Adjust perplexity according to the number of data points
        return TSNE(
            n_components=n_components,
            perplexity=perplexity,
            random_state=np.random.RandomState(seed),
        ).fit_transform(X)

    def do_bhsne(n_components=n_components, perplexity=perplexity, seed=42):
        return bh_sne(
            data=X,
            d=n_components,
            perplexity=perplexity,
            random_state=np.random.RandomState(seed),
        )

    def do_UMAP(n_neighbors=15, n_components=n_components, metric="euclidean"):
        return UMAP(
            n_neighbors=n_neighbors, n_components=n_components, metric=metric
        ).fit_transform(X)

    dispatcher = {"sksne": do_sksne, "bhsne": do_bhsne, "umap": do_UMAP}
    logger.debug(f"Performing embedding with {method}")
    X = dispatcher[method](**kwargs)
    if n_components == 3:
        embedded_df = pd.DataFrame(X, columns=["x", "y", "z"], index=df.index)
    elif n_components == 2:
        embedded_df = pd.DataFrame(X, columns=["x", "y"], index=df.index)
    else:
        embedded_df = pd.DataFrame(X, index=df.index)
    if embedded:
        embedded_df.to_csv(embedded, sep="\t", index=True, header=True)
        logger.debug(f"embedded.shape {embedded_df.shape} : Written {embedded}")
    return embedded_df


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    skip_desc = "(will skip if file exists)"
    cpus = mp.cpu_count()
    parser = argparse.ArgumentParser(
        description="Count k-mer frequencies of given `fasta`"
    )
    parser.add_argument(
        "--fasta", help="Metagenomic assembly fasta file", required=True
    )
    parser.add_argument(
        "--kmers",
        help=f"K-mers frequency tab-delimited table {skip_desc}",
        required=True,
    )
    parser.add_argument("--size", help="k-mer size in bp", default=5, type=int)
    parser.add_argument(
        "--normalized", help=f"</path/to/output/kmers.normalized.tsv> {skip_desc}"
    )
    parser.add_argument(
        "--embedded", help=f"</path/to/output/kmers.embedded.tsv> {skip_desc}"
    )
    parser.add_argument(
        "--method",
        help="embedding method [sk,bh]sne are corresponding implementations from scikit-learn and tsne, respectively.",
        choices=["sksne", "bhsne", "umap"],
        default="bhsne",
    )
    parser.add_argument(
        "--n-components",
        help="Number of dimensions to reduce k-mer frequencies to",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--do-pca",
        help="Whether to perform PCA prior to dimension reduction",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--pca-dimensions",
        help="<num components to reduce to PCA feature space",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--multiprocess",
        help="count k-mers using multiprocessing",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--nproc",
        help=f"num. processors to use if multiprocess is selected. (default = {cpus})",
        default=cpus,
        type=int,
    )
    args = parser.parse_args()
    try:
        df = load(args.kmers)
        logger.debug(f"{args.kmers} exists... loaded: df.shape {df.shape}")
    except FileNotFoundError as err:
        df = count(
            assembly=args.fasta,
            kmer_size=args.size,
            normalized=False,
            multiprocess=args.multiprocess,
            nproc=args.nproc,
        )
        df.to_csv(args.kmers, sep="\t", header=True, index=True)
        logger.debug(
            f"Wrote {len(df)} contigs {args.size}-mers frequencies to {args.kmers}."
        )

    if args.normalized:
        ndf = None
        try:
            ndf = load(args.normalized)
            logger.debug(f"{args.normalized} exists... loaded: df.shape {ndf.shape}")
        except FileNotFoundError as err:
            logger.debug(f"Normalizing {df.shape} k-mers DataFrame.")
            ndf = normalize(df)
            ndf.to_csv(args.normalized, sep="\t", header=True, index=True)
            logger.debug(
                f"Wrote {len(df)} normalized k-mer freqs. to {args.normalized}."
            )

    if not args.embedded:
        return

    if args.normalized:
        logger.debug(f"Embedding {args.normalized}")
    else:
        logger.debug(f"Embedding {args.kmers}")

    kmers_fp = args.normalized if args.normalized else args.kmers
    embedded_df = embed(
        kmers=kmers_fp,
        embedded=args.embedded,
        method=args.method,
        n_components=args.n_components,
        do_pca=args.do_pca,
        pca_dimensions=args.pca_dimensions,
    )


if __name__ == "__main__":
    main()
