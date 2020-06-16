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
Functions related to running hmmer on metagenome sequences
"""


import logging
import os
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

from glob import glob
from Bio import SeqIO

from autometa.common.external import prodigal


logger = logging.getLogger(__name__)


def annotate_parallel(orfs, hmmdb, outfpath, cpus):
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
    proc = subprocess.run(cmdline, shell=True)
    try:
        proc.check_returncode()
    except subprocess.CalledProcessError as err:
        logger.warning(f"Make sure your hmm profiles are pressed! hmmpress -f {hmmdb}")
        logger.error(f"Args:{cmdline} ReturnCode:{proc.returncode}")
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


def annotate_sequential(orfs, hmmdb, outfpath, cpus):
    cmd = ["hmmscan", "--cpu", str(cpus), "--tblout", outfpath, hmmdb, orfs]
    logger.debug(" ".join(cmd))
    proc = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        proc.check_returncode()
    except subprocess.CalledProcessError as err:
        logger.warning(f"Make sure your hmm profiles are pressed! hmmpress -f {hmmdb}")
        logger.error(f"Args:{cmd} ReturnCode:{proc.returncode}")
        raise err


def hmmscan(orfs, hmmdb, outfpath, cpus=0, force=False, parallel=True):
    """Runs hmmscan on dataset ORFs and provided hmm database.

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
        Will parallelize hmmscan using GNU parallel (the default is True).
    log : str, optional
        </path/to/parallel.log> (the default is None). If provided will write
        parallel log to `log`.

    Returns
    -------
    str
        </path/to/output.hmmscan.tsv>

    Raises
    -------
    FileExistsError
        `outfpath` already exists
    subprocess.CalledProcessError
        hmmscan failed
    """
    # OPTIMIZE: we want to extend parallel to grid computing (workqueue?) via --sshlogin?
    if os.path.exists(outfpath) and os.path.getsize(outfpath) > 0 and not force:
        raise FileExistsError(f"{outfpath}. Use force to overwrite!")
    if parallel:
        annotate_parallel(orfs=orfs, hmmdb=hmmdb, outfpath=outfpath, cpus=cpus)
    else:
        annotate_sequential(orfs=orfs, hmmdb=hmmdb, outfpath=outfpath, cpus=cpus)
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
    proc = subprocess.run(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )
    return fpath


def filter_markers(infpath, outfpath, cutoffs, orfs=None, force=False):
    """Filter markers from hmmscan output table that are above cutoff values.

    Parameters
    ----------
    infpath : str
        </path/to/hmmscan.tsv>
    outfpath : str
        </path/to/output.markers.tsv>
    cutoffs : str
        </path/to/cutoffs.tsv>
    orfs : str, optional
        Default will attempt to translate recovered qseqids to contigs
        </path/to/prodigal/called/orfs.fasta>
    force : bool, optional
        Overwrite existing `outfpath` (the default is False).

    Returns
    -------
    str
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
    if os.path.exists(outfpath) and os.path.getsize(outfpath) > 0 and not force:
        raise FileExistsError(f"{outfpath} already exists")
    hmmtab_header = ["sname", "sacc", "orf", "score"]
    col_indices = [0, 1, 2, 5]
    columns = {i: k for i, k in zip(col_indices, hmmtab_header)}
    df = pd.read_csv(infpath, sep="\s+", usecols=col_indices, header=None, comment="#")
    df.rename(columns=columns, inplace=True)
    # NaN may result from parsing issues while merging parallel results
    df.dropna(inplace=True)
    df["cleaned_sacc"] = df["sacc"].map(lambda acc: acc.split(".")[0])
    dff = pd.read_csv(cutoffs, sep="\t", index_col="accession")
    mdf = pd.merge(df, dff, how="left", left_on="cleaned_sacc", right_on="accession")
    mdf = mdf[mdf["score"] >= mdf["cutoff"]]
    if mdf.empty:
        raise AssertionError(f"No markers in {infpath} pass cutoff thresholds")
    cols = ["orf", "sacc", "sname", "score", "cutoff"]
    mdf = mdf[cols]
    translations = prodigal.contigs_from_headers(orfs)

    def translater(x):
        return translations.get(x, x.rsplit("_", 1)[0])

    mdf["contig"] = mdf["orf"].map(translater)
    mdf.set_index("contig", inplace=True)
    mdf.to_csv(outfpath, sep="\t", index=True, header=True)
    return outfpath


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
    parser.add_argument("hmmscan", help="</path/to/hmmscan.out>")
    parser.add_argument("markers", help="</path/to/markers.tsv>")
    parser.add_argument(
        "--force", help="force overwrite of out filepath", action="store_true"
    )
    parser.add_argument("--cpus", help="num cpus to use", default=0, type=int)
    parser.add_argument("--parallel", help="enable GNU parallel", action="store_true")
    args = parser.parse_args()

    if (
        os.path.exists(args.hmmscan)
        and os.path.getsize(args.hmmscan) > 0
        and not args.force
    ):
        result = args.hmmscan
    else:
        result = hmmscan(
            orfs=args.orfs,
            hmmdb=args.hmmdb,
            outfpath=args.hmmscan,
            cpus=args.cpus,
            force=args.force,
            parallel=args.parallel,
        )

    result = filter_markers(
        infpath=result,
        outfpath=args.markers,
        cutoffs=args.cutoffs,
        orfs=args.orfs,
        force=args.force,
    )


if __name__ == "__main__":
    main()
