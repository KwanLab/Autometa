#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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
Functions related to running hmmer on metagenome sequences
"""


import logging
import os
import shutil
import subprocess
import sys
import tempfile

import multiprocessing as mp
import pandas as pd

from glob import glob

from autometa.common.external import prodigal


logger = logging.getLogger(__name__)


def annotate_parallel(orfs, hmmdb, outfpath, cpus, seed=42):
    outdir = os.path.dirname(os.path.realpath(outfpath))
    outprefix = os.path.splitext(os.path.basename(outfpath))[0]
    tmp_dirpath = tempfile.mkdtemp(dir=outdir)
    __, tmp_fpath = tempfile.mkstemp(
        suffix=".{#}.txt", prefix=outprefix, dir=tmp_dirpath
    )
    log = os.path.join(outdir, "hmmscan.parallel.log")
    jobs = f"-j{cpus}"
    cmd = [
        "parallel",
        "--retries",
        "4",
        "--joblog",
        log,
        jobs,
        "--linebuffer",
        "--pipe",
        "--recstart",
        "'>'",
        "hmmscan",
        "--seed",
        str(seed),
        "--cpu",
        "0",
        "-o",
        os.devnull,
        "--tblout",
        tmp_fpath,
        hmmdb,
        "-",
        "<",
        orfs,
        "2>",
        os.devnull,
    ]
    cmdline = subprocess.list2cmdline(cmd)
    logger.debug(cmdline)
    try:
        subprocess.run(cmdline, shell=True, check=True)
    except subprocess.CalledProcessError as err:
        logger.warning(f"Make sure your hmm profiles are pressed! hmmpress -f {hmmdb}")
        shutil.rmtree(tmp_dirpath)
        raise err
    tmp_fpaths = glob(os.path.join(tmp_dirpath, "*.txt"))
    lines = ""
    buffer_limit = 60000
    with open(outfpath, "w") as out:
        for fp in tmp_fpaths:
            with open(fp) as fh:
                for line in fh:
                    lines += line
                    if sys.getsizeof(lines) >= buffer_limit:
                        out.write(lines)
                        lines = ""
        out.write(lines)
    shutil.rmtree(tmp_dirpath)


def annotate_sequential(orfs, hmmdb, outfpath, cpus, seed=42):
    cmd = [
        "hmmscan",
        "--seed",
        str(seed),
        "--cpu",
        str(cpus),
        "--tblout",
        outfpath,
        hmmdb,
        orfs,
    ]
    logger.debug(" ".join(cmd))
    try:
        subprocess.run(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
        )
    except subprocess.CalledProcessError as err:
        logger.warning(f"Make sure your hmm profiles are pressed! hmmpress -f {hmmdb}")
        raise err


def run(
    orfs,
    hmmdb,
    outfpath,
    cpus=0,
    force=False,
    parallel=True,
    gnu_parallel=False,
    seed=42,
):
    """Runs hmmscan on dataset ORFs and provided hmm database.

    Note
    ----
    Only one of `parallel` and `gnu_parallel` may be provided as True

    Parameters
    ----------
    orfs : str
        </path/to/orfs.faa>
    hmmdb : str
        </path/to/hmmpressed/database.hmm>
    outfpath : str
        </path/to/output.hmmscan.tsv>
    cpus : int, optional
        Num. cpus to use. 0 will run as many cpus as possible (the default is 0).
    force : bool, optional
        Overwrite existing `outfpath` (the default is False).
    parallel : bool, optional
        Will use multithreaded parallelization offered by hmmscan (the default is True).
    gnu_parallel : bool, optional
        Will parallelize hmmscan using GNU parallel (the default is False).
    seed : int, optional
        set RNG seed to <n> (if 0: one-time arbitrary seed) (the default is 42).

    Returns
    -------
    str
        </path/to/output.hmmscan.tsv>

    Raises
    -------
    ValueError
        Both `parallel` and `gnu_parallel` were provided as True
    FileExistsError
        `outfpath` already exists
    subprocess.CalledProcessError
        hmmscan failed

    """
    if gnu_parallel and parallel:
        raise ValueError("Both parallel and gnu_parallel were provided as True")
    # OPTIMIZE: we want to extend parallel to grid computing (workqueue?) via --sshlogin?
    if os.path.exists(outfpath) and os.path.getsize(outfpath) and not force:
        raise FileExistsError(f"{outfpath}. Use force to overwrite!")
    if gnu_parallel:
        annotate_parallel(
            orfs=orfs, hmmdb=hmmdb, outfpath=outfpath, cpus=cpus, seed=seed
        )
    elif parallel:
        cpus = mp.cpu_count() if not cpus else cpus
        annotate_sequential(
            orfs=orfs, hmmdb=hmmdb, outfpath=outfpath, cpus=cpus, seed=seed
        )
    else:
        annotate_sequential(
            orfs=orfs, hmmdb=hmmdb, outfpath=outfpath, cpus=0, seed=seed
        )
    if not os.path.exists(outfpath):
        raise FileNotFoundError(f"{outfpath} not written.")
    return outfpath


def hmmpress(fpath):
    """Runs hmmpress on `fpath`.

    Parameters
    ----------
    fpath : str
        </path/to/kindom.markers.hmm>

    Returns
    -------
    str
        </path/to/hmmpressed/kindom.markers.hmm>

    Raises
    -------
    FileNotFoundError
        `fpath` not found.
    subprocess.CalledProcessError
        hmmpress failed
    """
    if not os.path.exists(fpath):
        raise FileNotFoundError(fpath)
    cmd = ["hmmpress", "-f", fpath]
    logger.debug(" ".join(cmd))
    subprocess.run(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )
    return fpath


def read_tblout(infpath: str) -> pd.DataFrame:
    """Read hmmscan tblout into pd.DataFrame

    Parameters
    ----------
    infpath : str
        Path to hmmscan_results.tblout

    Returns
    -------
    pd.DataFrame
        DataFrame of raw hmmscan results

    Raises
    ------
    FileNotFoundError
        Path to `infpath` was not found
    """
    if not os.path.exists(infpath):
        raise FileNotFoundError(infpath)
    return pd.read_csv(
        infpath,
        sep="\s+",
        usecols=range(0, 18),
        names=[
            "sname",
            "sacc",
            "qname",
            "qacc",
            "full_seq_evalue",
            "full_seq_score",
            "full_seq_bias",
            "best_domain_evalue",
            "best_domain_score",
            "best_domain_bias",
            "domain_number_exp",
            "domain_number_reg",
            "domain_number_clu",
            "domain_number_ov",
            "domain_number_env",
            "domain_number_dom",
            "domain_number_rep",
            "domain_number_inc",
        ],
        comment="#",
    )


def read_domtblout(fpath: str) -> pd.DataFrame:
    """Read hmmscan domtblout file into pandas DataFrame

    For more detailed column descriptions see the 'tabular output formats' section in the [HMMER manual](http://eddylab.org/software/hmmer/Userguide.pdf#tabular-output-formats)

    Parameters
    ----------
    fpath : str
        Path to hmmscan domtblout file

    Returns
    -------
    pd.DataFrame
        index=range(0,n_hits) cols=...
    """
    if not os.path.exists(fpath):
        raise FileNotFoundError(fpath)
    return pd.read_csv(
        fpath,
        sep="\s+",
        usecols=range(0, 22),
        names=[
            "sname",
            "sacc",
            "slen",
            "qname",
            "qacc",
            "qlen",
            "full_seq_evalue",
            "full_seq_score",
            "full_seq_bias",
            "this_domain_count",
            "total_domain_count",
            "domain_conditional_evalue",
            "domain_independent_evalue",
            "domain_score",
            "domain_bias",
            "hmm_coord_from",
            "hmm_coord_to",
            "ali_coord_from",
            "ali_coord_to",
            "env_coord_from",
            "env_coord_to",
            "mean_posterior_probability_of_aligned_residues",
        ],
        comment="#",
    )


def filter_markers(
    infpath: str,
    cutoffs: str,
    outfpath: str = None,
    orfs: str = None,
    force: bool = False,
) -> pd.DataFrame:
    """Filter markers from hmmscan tblout output table using provided cutoff values file.

    Parameters
    ----------
    infpath : str
        Path to hmmscan tblout output file
    cutoffs : str
        Path to marker set inclusion cutoffs
    outfpath : str, optional
        Path to write filtered markers to tab-delimited file
    orfs : str, optional
        Default will attempt to translate recovered qseqids to contigs
        </path/to/prodigal/called/orfs.fasta>
    force : bool, optional
        Overwrite existing `outfpath` (the default is False).

    Returns
    -------
    pd.DataFrame
        </path/to/output.markers.tsv>

    Raises
    -------
    FileNotFoundError
        `infpath` or `cutoffs` not found
    FileExistsError
        `outfpath` already exists and force=False
    AssertionError
        No returned markers pass the cutoff thresholds. I.e. final df is empty.
    """
    for fp in [infpath, cutoffs]:
        if not os.path.exists(fp):
            raise FileNotFoundError(fp)
    if os.path.exists(outfpath) and os.path.getsize(outfpath) and not force:
        raise FileExistsError(f"{outfpath} already exists")

    df = read_tblout(infpath).dropna()
    df["cleaned_sacc"] = df.sacc.map(lambda acc: acc.split(".")[0])
    dff = pd.read_csv(cutoffs, sep="\t", index_col="accession")
    mdf = pd.merge(df, dff, how="left", left_on="cleaned_sacc", right_on="accession")
    # Filter markers using cutoffs compared to full seq score
    mdf = mdf.loc[mdf.full_seq_score.ge(mdf.cutoff)]
    if mdf.empty:
        raise AssertionError(f"No markers in {infpath} pass cutoff thresholds")

    translations = prodigal.contigs_from_headers(orfs) if orfs else {}

    def translater(x):
        return translations.get(x, x.rsplit("_", 1)[0])

    mdf["contig"] = mdf.qname.map(translater)
    mdf.set_index("contig", inplace=True)
    if outfpath:
        mdf.to_csv(outfpath, sep="\t", index=True, header=True)
    return mdf


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Retrieves markers with provided input assembly"
    )
    parser.add_argument("orfs", help="</path/to/assembly.orfs.faa>")
    parser.add_argument("hmmdb", help="</path/to/hmmpressed/hmmdb>")
    parser.add_argument("cutoffs", help="</path/to/hmm/cutoffs.tsv>")
    parser.add_argument("hmmscan", help="</path/to/hmmscan.tblout>")
    parser.add_argument("markers", help="</path/to/markers.tsv>")
    parser.add_argument(
        "--force", help="force overwrite of out filepath", action="store_true"
    )
    parser.add_argument("--cpus", help="num cpus to use", default=0, type=int)
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--parallel",
        help="enable hmmer multithreaded parallelization",
        action="store_true",
    )
    group.add_argument(
        "--gnu-parallel", help="enable GNU parallelization", action="store_true"
    )
    args = parser.parse_args()

    if (
        os.path.exists(args.hmmscan)
        and os.path.getsize(args.hmmscan)
        and not args.force
    ):
        result = args.hmmscan
    else:
        result = run(
            orfs=args.orfs,
            hmmdb=args.hmmdb,
            outfpath=args.hmmscan,
            cpus=args.cpus,
            force=args.force,
            parallel=args.parallel,
            gnu_parallel=args.gnu_parallel,
        )

    df = filter_markers(
        infpath=result,
        outfpath=args.markers,
        cutoffs=args.cutoffs,
        orfs=args.orfs,
        force=args.force,
    )

    markers_df = df.groupby("contig")["sacc"].value_counts().unstack().convert_dtypes()


if __name__ == "__main__":
    main()
