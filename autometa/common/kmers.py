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
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
import multiprocessing as mp

from tqdm import tqdm
from Bio import SeqIO
from scipy.stats import gmean
from skbio.stats.composition import ilr, clr, multiplicative_replacement
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from tsne import bh_sne
from umap import UMAP

from autometa.common import utilities
from autometa.common.exceptions import TableFormatError

# Suppress numba logger debug output
numba_logger = logging.getLogger("numba").setLevel(logging.ERROR)
logger = logging.getLogger(__name__)


def _revcomp(string: str) -> str:
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


def init_kmers(kmer_size: int = 5) -> Dict[str, int]:
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


def load(kmers_fpath: str) -> pd.DataFrame:
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
    TableFormatError
        `kmers_fpath` file format is invalid

    """
    if not os.path.exists(kmers_fpath) or not os.path.getsize(kmers_fpath):
        raise FileNotFoundError(kmers_fpath)
    try:
        df = pd.read_csv(kmers_fpath, sep="\t", index_col="contig")
    except ValueError:
        raise TableFormatError(f"contig column not found in {kmers_fpath}!")
    return df


def mp_counter(
    assembly: str, ref_kmers: Dict[str, int], cpus: int = mp.cpu_count()
) -> List:
    """Multiprocessing k-mer counter used in `count`. (Should not be used directly).

    Parameters
    ----------
    assembly : str
        </path/to/assembly.fasta> (nucleotides)
    ref_kmers : dict
        {kmer:index, ...}
    cpus : int, optional
        Number of cpus to use. (the default will use all available).

    Returns
    -------
    list
        [{record:counts}, {record:counts}, ...]

    """
    pool = mp.Pool(cpus)
    args = [(record, ref_kmers) for record in SeqIO.parse(assembly, "fasta")]
    logger.debug(
        f"Pool counter (cpus={cpus}): counting {len(args):,} records k-mer frequencies"
    )
    results = pool.map(record_counter, args)
    pool.close()
    pool.join()
    counts = {}
    for result in results:
        counts.update(result)
    return counts


def record_counter(
    args: Tuple[SeqIO.SeqRecord, Dict[str, int]]
) -> Dict[str, List[int]]:
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


def seq_counter(
    assembly: str, ref_kmers: Dict[str, int], verbose: bool = True
) -> Dict[str, List[int]]:
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


@utilities.timeit
def count(
    assembly: str,
    size: int = 5,
    out: str = None,
    force: bool = False,
    verbose: bool = True,
    multiprocess: bool = True,
    cpus: int = mp.cpu_count(),
) -> pd.DataFrame:
    """Counts k-mer frequencies for provided assembly file

    First we make a dictionary of all the possible k-mers (discounting reverse
    complements). Each k-mer's count is updated by index when encountered in the
    record.

    Parameters
    ----------
    assembly : str
        Description of parameter `assembly`.
    size : int, optional
        length of k-mer to count `size` (the default is 5).
    out: str, optional
        Path to write k-mer counts table.
    force: bool, optional
        Whether to overwrite existing `out` k-mer counts table (the default is False).
    verbose : bool, optional
        Enable progress bar `verbose` (the default is True).
    multiprocess : bool, optional
        Use multiple cores to count k-mer frequencies (the default is True).
    cpus : int, optional
        Number of cpus to use. (the default will use all available).

    Returns
    -------
    pandas.DataFrames
        index_col='contig', tab-delimited, cols=unique_kmers
        i.e. 5-mer columns=[AAAAA, AAAAT, AAAAC, AAAAG,  ..., GCCGC]

    Raises
    -------
    TypeError
        `size` must be an int

    """
    out_specified = out is not None
    out_exists = os.path.exists(out) if out else False
    case1 = out_specified and out_exists and not force
    if case1:
        logger.warning(f"counts already exist: {out} force to overwrite. [retrieving]")
        df = pd.read_csv(out, sep="\t", index_col="contig")
    else:
        if not isinstance(size, int):
            raise TypeError(f"size must be an int! Given: {type(size)}")
        ref_kmers = init_kmers(size)
        if assembly.endswith(".gz"):
            assembly = utilities.gunzip(assembly, assembly.rstrip(".gz"))
        logger.info(f"Counting {size}-mers.")
        if multiprocess:
            kmer_counts = mp_counter(assembly, ref_kmers, cpus)
        else:
            kmer_counts = seq_counter(assembly, ref_kmers, verbose=verbose)
        df = pd.DataFrame(kmer_counts, index=ref_kmers).transpose()
        df.index.name = "contig"
    if out:
        df.to_csv(out, sep="\t", index=True, header=True)
        logger.debug(f"Wrote {df.shape[0]} contigs {size}-mers frequencies to {out}.")
    return df


def autometa_clr(df: pd.DataFrame) -> pd.DataFrame:
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


def normalize(
    df: pd.DataFrame, method: str = "am_clr", out: str = None, force: bool = False
) -> pd.DataFrame:
    """Normalize raw k-mer counts by center or isometric log-ratio transform.

    Parameters
    ----------
    df : pd.DataFrame
        k-mer counts dataframe.
        i.e. for 3-mers; Index='contig', columns=[AAA, AAT, ...]
    method : str, optional
        Normalize k-mer counts using CLR or ILR transformation
        (the default is Autometa's CLR implementation).
        choices = ['ilr', 'clr', 'am_clr']
        Other transformations come from the skbio.stats.composition module
    out : str, optional
        Path to write normalized k-mers.
    force : bool, optional
        Whether to overwrite existing `out` file path, by default False.

    Returns
    -------
    pd.DataFrame
        Normalized counts using provided `method`.

    Raises
    ------
    ValueError
        Provided `method` is not available.
    """
    method = method.lower()
    out_specified = out is not None
    out_exists = os.path.exists(out) if out else False
    case1 = out_specified and out_exists and not force
    if case1:
        logger.debug(f"{out} already exists. Use force to overwrite. retrieving...")
        return pd.read_csv(out, sep="\t", index_col="contig")
    logger.debug(f"Transforming k-mer counts using {method}")
    choices = {"ilr", "clr", "am_clr"}
    if method == "am_clr":
        norm_df = autometa_clr(df)
    elif method in choices:
        transforms = {"ilr": ilr, "clr": clr}
        X = df.fillna(0).to_numpy()
        X = multiplicative_replacement(X)
        X_norm = transforms[method](X)
        norm_df = pd.DataFrame(X_norm, index=df.index)
    else:
        choices = ", ".join(choices)
        raise ValueError(
            f"Normalize Method not available! {method}. choices: {choices}"
        )
    case2 = out_specified and out_exists and force
    case3 = out_specified and not out_exists
    if case2 or case3:
        norm_df.to_csv(out, sep="\t", index=True, header=True)
    return norm_df


def embed(
    kmers: Union[str, pd.DataFrame],
    out: str = None,
    force: bool = False,
    embed_dimensions: int = 2,
    do_pca: bool = True,
    pca_dimensions: int = 50,
    method: str = "bhsne",
    perplexity: float = 30.0,
    seed: int = 42,
    **method_args: Dict,
) -> pd.DataFrame:
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
    out : str, optional
        </path/to/output/kmers.out.tsv> If provided will write to `out`.
    force: bool, optional
        Whether to overwrite existing `out` file.
    embed_dimensions : int, optional
        embedn_dimensions` to embed k-mer frequencies (the default is 2).
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
    seed: int, optional
        Seed to use for `method`. Allows for reproducibility from random state.
    **method_args : dict, optional
        Other arguments to be supplied to respective `method`.

    Returns
    -------
    pd.DataFrame
        out dataframe with index='contig' and cols=['x','y','z']

    Raises
    -------
    TypeError
        Provided `kmers` is not a str or pd.DataFrame.
    TableFormatError
        Provided `kmers` or `out` are not formatted correctly for use.
    ValueError
        Provided `method` is not an available choice.
    FileNotFoundError
        `kmers` type must be a pd.DataFrame or filepath.
    """
    if isinstance(kmers, str) and os.path.exists(kmers) and os.path.getsize(kmers):
        try:
            df = pd.read_csv(kmers, sep="\t", index_col="contig")
        except ValueError:
            raise TableFormatError(f"contig column not found in {kmers}")
    elif isinstance(kmers, pd.DataFrame):
        df = kmers
    else:
        raise TypeError(kmers)
    if out and os.path.exists(out) and os.path.getsize(out) and not force:
        logger.debug(f"k-mers frequency embedding already exists {out}")
        try:
            return pd.read_csv(out, sep="\t", index_col="contig")
        except ValueError:
            raise TableFormatError(f"contig column not found in {out}")

    if df.empty:
        kmers_desc = f"kmers:{kmers} type:{type(kmers)}"
        embed_desc = f"out:{out} type:{type(out)}"
        requirements = f"kmers type must be a pd.DataFrame or filepath."
        raise FileNotFoundError(f"{kmers_desc} {embed_desc} {requirements}")

    method = method.lower()
    choices = {"umap", "sksne", "bhsne"}
    if method not in choices:
        raise ValueError(
            f"{method} not in embedding methods. Choices: {', '.join(choices)}"
        )
    # PCA
    n_samples, n_components = df.shape
    # Drop any rows that all cols contain NaN. This may occur if the contig length is below the k-mer size
    df.dropna(axis="index", how="all", inplace=True)
    df.fillna(0, inplace=True)
    X = df.to_numpy()
    # Set random state using provided seed
    random_state = np.random.RandomState(seed)

    if n_components > pca_dimensions and do_pca:
        logger.debug(
            f"Performing decomposition with PCA (seed {seed}): {n_components} to {pca_dimensions} dims"
        )
        X = PCA(n_components=pca_dimensions, random_state=random_state).fit_transform(X)
        # X = PCA(n_components='mle').fit_transform(X)
        n_samples, n_components = X.shape

    logger.debug(f"{method}: {n_samples} data points and {n_components} dimensions")

    # Adjust perplexity according to the number of data points
    n_rows = n_samples - 1
    scaler = 3.0
    if n_rows < (scaler * perplexity):
        perplexity = (n_rows / scaler) - 1

    def do_sksne():
        return TSNE(
            n_components=embed_dimensions,
            perplexity=perplexity,
            random_state=random_state,
        ).fit_transform(X)

    def do_bhsne():
        return bh_sne(
            data=X, d=embed_dimensions, perplexity=perplexity, random_state=random_state
        )

    def do_UMAP():
        return UMAP(
            n_neighbors=15,
            n_components=embed_dimensions,
            metric="euclidean",
            random_state=random_state,
        ).fit_transform(X)

    dispatcher = {"sksne": do_sksne, "bhsne": do_bhsne, "umap": do_UMAP}
    logger.debug(f"Performing embedding with {method} (seed {seed})")
    try:
        X = dispatcher[method](**method_args)
    except ValueError as err:
        if method == "sksne":
            logger.warning(
                f"--embed-dimensions ({embed_dimensions}) is too high for sksne. Reducing to 3."
            )
            embed_dimensions = 3
            X = dispatcher[method](**method_args)
        else:
            raise err
    if embed_dimensions == 3:
        embedded_df = pd.DataFrame(X, columns=["x", "y", "z"], index=df.index)
    elif embed_dimensions == 2:
        embedded_df = pd.DataFrame(X, columns=["x", "y"], index=df.index)
    else:
        embedded_df = pd.DataFrame(X, index=df.index)
    if out:
        embedded_df.to_csv(out, sep="\t", index=True, header=True)
        logger.debug(f"embedded.shape {embedded_df.shape} : Written {out}")
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
        "--norm-method",
        help="embedding method [sk,bh]sne are corresponding implementations from scikit-learn and tsne, respectively.",
        choices=["ilr", "clr", "am_clr"],
        default="am_clr",
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
        "--embedded", help=f"</path/to/output/kmers.embedded.tsv> {skip_desc}"
    )
    parser.add_argument(
        "--embed-method",
        help="embedding method [sk,bh]sne are corresponding implementations from scikit-learn and tsne, respectively.",
        choices=["sksne", "bhsne", "umap"],
        default="bhsne",
    )
    parser.add_argument(
        "--embed-dimensions",
        help="Number of dimensions to reduce k-mer frequencies to",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--force",
        help="Whether to overwrite existing annotations",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--multiprocess",
        help="count k-mers using multiprocessing",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--cpus",
        help=f"num. processors to use if multiprocess is selected. (default = {cpus})",
        default=cpus,
        type=int,
    )
    parser.add_argument(
        "--seed",
        help=f"Seed to set random state for dimension reduction determinism.",
        default=42,
        type=int,
    )
    args = parser.parse_args()

    df = count(
        assembly=args.fasta,
        size=args.size,
        out=args.kmers,
        force=args.force,
        multiprocess=args.multiprocess,
        cpus=args.cpus,
    )

    if args.normalized:
        df = normalize(
            df=df, method=args.norm_method, out=args.normalized, force=args.force
        )

    if args.embedded:
        embedded_df = embed(
            kmers=df,
            out=args.embedded,
            force=args.force,
            method=args.embed_method,
            embed_dimensions=args.embed_dimensions,
            do_pca=args.do_pca,
            pca_dimensions=args.pca_dimensions,
            seed=args.seed,
        )


if __name__ == "__main__":
    main()
