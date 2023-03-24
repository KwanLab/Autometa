#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Count, normalize and embed k-mers given nucleotide sequences
"""


import gzip
import logging
import os
import sys
from typing import Any, Dict, List, Tuple, Union

import numpy as np
import pandas as pd
import multiprocessing as mp

from tqdm import tqdm
from Bio import SeqIO
from scipy.stats import gmean
from autometa.common.file_handling import open_file
from skbio.stats.composition import ilr, clr, multiplicative_replacement
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from tsne import bh_sne
from umap import UMAP
from trimap import TRIMAP

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
    if not isinstance(kmer_size, int):
        raise TypeError(f"kmer_size must be an int! Given: {type(kmer_size)}")
    index = 0
    uniq_kmers = dict()
    dna_letters = ["A", "T", "C", "G"]
    all_kmers = dna_letters
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
        Path to kmer frequency table

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
    except (ValueError, TypeError):
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
    if assembly.endswith(".gz"):
        fh = gzip.open(assembly, "rt")
        args = [(record, ref_kmers) for record in SeqIO.parse(fh, "fasta")]
        fh.close()
    else:
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
    fh = gzip.open(assembly, "rt") if assembly.endswith(".gz") else assembly
    records = SeqIO.parse(fh, "fasta")
    desc = f"Counting {n_uniq_kmers} unique {kmer_size}-mers"
    disable = not verbose
    for record in tqdm(records, desc=desc, disable=disable, leave=False):
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
    if hasattr(fh, "close"):
        # Ensure we close the open gzipped file
        fh.close()
    return kmer_counts


@utilities.timeit
def count(
    assembly: str,
    size: int = 5,
    out: str = None,
    force: bool = False,
    verbose: bool = True,
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
    if out_specified and out_exists and not force:
        logger.warning(f"counts already exist: {out} force to overwrite. [retrieving]")
        df = pd.read_csv(out, sep="\t", index_col="contig")
    else:
        # checks if it is something like 12.0 vs. 12.9. Also check is an int
        if not isinstance(size, int):
            raise TypeError(f"size must be an int! Given: {type(size)}")
        ref_kmers = init_kmers(size)
        logger.info(f"Counting {size}-mers.")
        if cpus > 1:
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

    Steps
    -----

        * Drop any k-mers not present for all contigs
        * Drop any contigs not containing any kmer counts
        * Fill any remaining na values with 0
        * Normalize the k-mer count by the total count of all k-mers for a given contig
        * Add 1 as 0 can not be utilized for CLR
        * Perform CLR transformation log(norm. value / geometric mean norm. value)

    Parameters
    ----------
    df : pd.DataFrame
        K-mers Dataframe where index_col='contig' and column values are k-mer
        frequencies.

    References
    ----------

        * Aitchison, J. The Statistical Analysis of Compositional Data (1986)
        * Pawlowsky-Glahn, Egozcue, Tolosana-Delgado. Lecture Notes on Compositional Data Analysis (2011)
        * Why ILR is preferred `stats stackexchange discussion <https://stats.stackexchange.com/questions/242445/why-is-isometric-log-ratio-transformation-preferred-over-the-additivealr-or-ce>`_
        * Use of CLR transformation prior to PCA `stats stackexchange discussion <https://stats.stackexchange.com/questions/305965/can-i-use-the-clr-centered-log-ratio-transformation-to-prepare-data-for-pca>`_
        * Lecture notes on Compositional Data Analysis (CoDa) `PDF <http://www.sediment.uni-goettingen.de/staff/tolosana/extra/CoDa.pdf>`_

    Returns
    -------
    pd.DataFrame
        index='contig', cols=[kmer, kmer, ...]
        Columns have been transformed by CLR normalization.

    """
    # steps in 1: data cleaning
    df = df.dropna(axis="columns", how="all").dropna(axis="index", how="all").fillna(0)
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
    pca_dimensions: int = 50,
    method: str = "bhsne",
    perplexity: float = 30.0,
    seed: int = 42,
    n_jobs: int = -1,
    **method_kwargs: Dict[str, Any],
) -> pd.DataFrame:
    """Embed k-mers using provided `method`.

    Notes
    -----

        * `sklearn.manifold.TSNE <https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html#sklearn.manifold.TSNE>`_
        * `tsne.bh_sne <https://pypi.org/project/tsne/>`_
        * `UMAP <https://umap-learn.readthedocs.io/en/latest/>`_
        * `densMAP <https://umap-learn.readthedocs.io/en/latest/densmap_demo.html#better-preserving-local-density-with-densmap>`_
        * `TriMap <https://github.com/eamid/trimap>`_

    Parameters
    ----------
    kmers : str or pd.DataFrame
        </path/to/input/kmers.normalized.tsv>

    out : str, optional
        </path/to/output/kmers.out.tsv> If provided will write to `out`.

    force: bool, optional
        Whether to overwrite existing `out` file.

    embed_dimensions : int, optional
        embed_dimensions` to embed k-mer frequencies (the default is 2).

        The output embedded kmers will follow columns of `x_1` to `x_{embed_dimensions}`

        NOTE: The columns are 1-indexed, i.e. at x_1 *not* x_0

    pca_dimensions : int, optional
        Reduce k-mer frequencies dimensions to `pca_dimensions` (the default is 50).
        If zero, will skip this step.

    method : str, optional
        embedding method to use (the default is 'bhsne').
        choices include sksne, bhsne, umap, trimap and densmap.

    perplexity : float, optional
        hyperparameter used to tune sksne and bhsne (the default is 30.0).

    seed: int, optional
        Seed to use for `method`. Allows for reproducibility from random state.

    n_jobs: int, optional

        Used with `sksne`, `densmap` and `umap`, (the default is -1 which will attempt to use all available CPUs)

        Note
        ----

        For n_jobs below -1, (CPUS + 1 + n_jobs) are used. For example with n_jobs=-2, all CPUs but one are used.

        * scikit-learn TSNE `n_jobs glossary <https://scikit-learn.org/stable/glossary.html#term-n_jobs>`_
        * UMAP and DensMAP's
        `invocation <https://github.com/lmcinnes/umap/blob/2c5232f7b946efab30e279c0b095b37f5648ed8b/umap/umap_.py#L328-L341>`_
        use this with
        `pynndescent <https://github.com/lmcinnes/pynndescent/blob/cc6ed32e25f7afb14913bff04d3b01723b33e5b5/pynndescent/pynndescent_.py#L629-L632>`_


    **method_kwargs : Dict[str, Any], optional

        Other keyword arguments (kwargs) to be supplied to respective `method`.

        Examples
        --------

        Set UMAP(verbose=True, output_dens=True) using **method_kwargs
        >>> embed_df = kmers.embed(
            norm_df,
            method='densmap',
            embed_dimensions=2,
            n_jobs=None,
            **{
                'verbose': True,
                'output_dens': True,
            }
        )

        NOTE: Setting duplicate arguments will result in an error

        Here we specify ``UMAP(densmap=True)`` using ``method='densmap'``
        and also attempt to overwrite to ``UMAP(densmap=False)``
        with the method_kwargs, ``**{'densmap':False}``, resulting
        in a TypeError.

        >>> embed_df = kmers.embed(
            df,
            method='densmap',
            embed_dimensions=2,
            n_jobs=4,
            **{'densmap': False}
        )
        TypeError: umap.umap_.UMAP() got multiple values for keyword argument 'densmap'

        Typically, you will not require the use of method_kwargs as this is only available
        for applying advanced parameter settings to any of the available embedding methods.

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
            with open_file(kmers) as h:
                df = pd.read_csv(h, sep="\t", index_col="contig")
        except ValueError:
            raise TableFormatError(f"contig column not found in {kmers}")
    elif isinstance(kmers, pd.DataFrame):
        df = kmers
    else:
        raise TypeError(kmers)
    if out and os.path.exists(out) and os.path.getsize(out) and not force:
        logger.debug(f"k-mers frequency embedding already exists {out}")
        try:
            with open_file(kmers) as h:
                return pd.read_csv(kmers, sep="\t", index_col="contig")
        except ValueError:
            raise TableFormatError(f"contig column not found in {out}")

    if df.empty:
        kmers_desc = f"kmers:{kmers} type:{type(kmers)}"
        embed_desc = f"out:{out} type:{type(out)}"
        requirements = f"Given pd.DataFrame is empty!"
        raise FileNotFoundError(f"{kmers_desc} {embed_desc} {requirements}")

    method = method.lower()
    choices = {"umap", "sksne", "bhsne", "densmap", "trimap"}
    if method not in choices:
        raise ValueError(
            f"{method} not in embedding methods. Choices: {', '.join(choices)}"
        )
    # PCA
    n_samples, n_components = df.shape
    # Drop any rows that all cols contain NaN. This may occur if the contig length is below the k-mer size
    X = df.dropna(axis="index", how="all").fillna(0).to_numpy()
    # Set random state using provided seed
    random_state = np.random.RandomState(seed)
    if isinstance(pca_dimensions, str):
        try:
            int(pca_dimensions)
        except Exception as e:
            raise TypeError(
                f"pca_dimensions must be an integer! given: {pca_dimensions}"
            )
    if df.shape[0] < pca_dimensions:
        logger.info(
            f"Stopping. Number of contigs ({str(df.shape[0])}) is less than pca_dimensions ({str(pca_dimensions)})"
        )
        # exit with 0 if not enough contigs are present
        # don't want to raise an actual error, because even if ignored Nextflow users will see a Note and be confused
        # TODO: write a file with the log message
        sys.exit(0)
    if n_components > pca_dimensions and pca_dimensions != 0:
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
            n_jobs=n_jobs,
            **method_kwargs,
        ).fit_transform(X)

    def do_bhsne():
        return bh_sne(
            data=X,
            d=embed_dimensions,
            perplexity=perplexity,
            random_state=random_state,
            **method_kwargs,
        )

    # def do_densne():
    #     return densne.run_densne(
    #         X, no_dims=embed_dimensions, perplexity=perplexity, rand_seed=random_state, **method_kwargs
    #     )

    method_is_densmap = method == "densmap"

    def do_UMAP():
        return UMAP(
            n_neighbors=15,
            n_components=embed_dimensions,
            metric="euclidean",
            random_state=random_state,
            densmap=method_is_densmap,
            n_jobs=n_jobs,
            **method_kwargs,
        ).fit_transform(X)

    def do_trimap():
        return TRIMAP(
            n_dims=embed_dimensions, verbose=False, **method_kwargs
        ).fit_transform(X)

    # TODO: Add "densne":do_densne() to dispatcher when easy install of densne is available.
    dispatcher = {
        "sksne": do_sksne,
        "bhsne": do_bhsne,
        "umap": do_UMAP,
        "densmap": do_UMAP,
        "trimap": do_trimap,
    }
    logger.debug(f"Performing embedding with {method} (seed {seed})")
    try:
        X = dispatcher[method]()
    except ValueError as err:
        if method == "sksne":
            logger.warning(
                f"embed_dimensions ({embed_dimensions}) is too high for sksne. Reducing to 3."
            )
            embed_dimensions = 3
            X = dispatcher[method]()
        else:
            raise err

    embed_cols = [f"x_{col}" for col in range(1, embed_dimensions + 1)]
    if isinstance(X, tuple):
        # When method_kwargs = **{'output_dens': True}
        # X : tuple[np.ndarray, np.ndarray, np.ndarray]
        # X : tuple[embedding, original local radii, embedding local radii]
        output_dens_ndarray_cols = [
            embed_cols,
            ["original_local_radius"],
            ["embedded_local_radius"],
        ]
        embedded_df = pd.concat(
            [
                pd.DataFrame(result, index=df.index, columns=cols)
                for result, cols in zip(X, output_dens_ndarray_cols)
            ],
            axis=1,
        )
    elif isinstance(X, np.ndarray):
        embedded_df = pd.DataFrame(X, index=df.index, columns=embed_cols)
    else:
        logger.warning(
            f"Unrecognized {method} transform (method_kwargs={method_kwargs}) output type: {type(X)}"
        )
        embedded_df = pd.DataFrame(X, index=df.index, columns=embed_cols)
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
    cpus = mp.cpu_count()
    parser = argparse.ArgumentParser(
        description="Count k-mer frequencies of given `fasta`",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--fasta",
        help="Metagenomic assembly fasta file",
        metavar="filepath",
    )
    parser.add_argument(
        "--kmers",
        help=f"K-mers frequency tab-delimited table (will skip if file exists)",
        metavar="filepath",
    )
    parser.add_argument(
        "--size", help="k-mer size in bp", default=5, metavar="int", type=int
    )
    parser.add_argument(
        "--norm-output",
        help=f"Path to normalized kmers table (will skip if file exists)",
        metavar="filepath",
    )
    parser.add_argument(
        "--norm-method",
        help="""Normalization method to transform kmer counts prior to PCA and embedding.
        ilr: isometric log-ratio transform (scikit-bio implementation).
        clr: center log-ratio transform (scikit-bio implementation).
        am_clr: center log-ratio transform (Autometa implementation).
        """,
        choices=["ilr", "clr", "am_clr"],
        default="am_clr",
    )
    parser.add_argument(
        "--pca-dimensions",
        help="Number of dimensions to reduce to PCA feature space after normalization and prior to embedding (NOTE: Setting to zero will skip PCA step)",
        type=int,
        metavar="int",
        default=50,
    )
    parser.add_argument(
        "--embedding-output",
        help=f"Path to write embedded kmers table (will skip if file exists)",
        metavar="filepath",
    )
    parser.add_argument(
        "--embedding-method",
        help="embedding method [sk,bh]sne are corresponding implementations from scikit-learn and tsne, respectively.",
        choices=["sksne", "bhsne", "umap", "densmap", "trimap"],
        default="bhsne",
    )
    parser.add_argument(
        "--embedding-dimensions",
        help="Number of dimensions of which to reduce k-mer frequencies",
        type=int,
        metavar="int",
        default=2,
    )
    parser.add_argument(
        "--force",
        help="Whether to overwrite existing annotations",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--cpus",
        help=f"num. processors to use.",
        default=cpus,
        metavar="int",
        type=int,
    )
    parser.add_argument(
        "--seed",
        help=f"Seed to set random state for dimension reduction determinism.",
        default=42,
        metavar="int",
        type=int,
    )
    args = parser.parse_args()

    if not args.fasta and not args.kmers and not args.norm_output:
        raise ValueError(
            "At least one of --fasta, --kmers or --norm-output are required!"
        )

    norm_df = pd.DataFrame()

    if (
        args.norm_output
        and not os.path.exists(args.norm_output)
        and not args.fasta
        and not args.kmers
    ):
        # only normalized kmers were provided
        raise FileNotFoundError(args.norm_output)
    elif args.kmers and not os.path.exists(args.kmers) and not args.fasta:
        # only kmer counts were provided
        raise FileNotFoundError(args.kmers)
    elif args.norm_output and os.path.exists(args.norm_output) and not args.force:
        # We already have the normalized kmers
        with open_file(args.norm_output) as h:
            norm_df = pd.read_csv(h, sep="\t", index_col="contig")
    elif args.kmers and os.path.exists(args.kmers) and not args.force:
        # We already have the kmer counts
        with open_file(args.kmers) as h:
            kmers_df = pd.read_csv(h, sep="\t", index_col="contig")
    else:
        # Start with counting kmers
        kmers_df = count(
            assembly=args.fasta,
            size=args.size,
            out=args.kmers,
            force=args.force,
            cpus=args.cpus,
        )

    if args.norm_output and norm_df.empty:
        norm_df = normalize(
            df=kmers_df,
            method=args.norm_method,
            out=args.norm_output,
            force=args.force,
        )

    if args.embedding_output:
        df = kmers_df if norm_df.empty else norm_df
        embedded_df = embed(
            kmers=df,
            out=args.embedding_output,
            force=args.force,
            method=args.embedding_method,
            embed_dimensions=args.embedding_dimensions,
            pca_dimensions=args.pca_dimensions,
            seed=args.seed,
            n_jobs=args.cpus,
        )


if __name__ == "__main__":
    main()
