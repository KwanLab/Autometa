#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Functions to retrieve orfs from provided assembly using prodigal
"""


import gzip
import logging
import os
import subprocess
import shutil
import tempfile

from glob import glob
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from typing import List, Set, Tuple, Union, Mapping

from autometa.config.environ import get_versions
from autometa.common.utilities import gunzip

logger = logging.getLogger(__name__)


def aggregate_orfs(search_str: str, outfpath: str) -> None:
    tmpfpaths = glob(search_str)
    lines = ""
    for fp in tmpfpaths:
        with open(fp) as fh:
            for line in fh:
                lines += line
    out = open(outfpath, "w")
    out.write(lines)
    out.close()


def annotate_sequential(assembly: str, prots_out: str, nucls_out: str) -> None:
    cmd = [
        "prodigal",
        "-i",
        assembly,
        "-a",
        prots_out,
        "-d",
        nucls_out,
        "-p",
        "meta",
        "-m",
        "-q",
    ]
    cmd = [str(arg) for arg in cmd]
    logger.debug(" ".join(cmd))
    subprocess.run(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )


def annotate_parallel(assembly: str, prots_out: str, nucls_out: str, cpus: int) -> None:
    outdir = os.path.dirname(os.path.realpath(nucls_out))
    log = os.path.join(outdir, "prodigal.parallel.log")
    outprefix = os.path.splitext(os.path.basename(nucls_out))[0]
    tmpdir = tempfile.mkdtemp(dir=outdir)
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)
    tmpnucl = ".".join([outprefix, "{#}", "fna"])
    tmpprot = ".".join([outprefix, "{#}", "faa"])
    tmpnucl_fpath = os.path.join(tmpdir, tmpnucl)
    tmpprot_fpath = os.path.join(tmpdir, tmpprot)
    jobs = f"-j{cpus}"
    cmd = [
        "parallel",
        "--retries",
        "4",
        "--joblog",
        log,
        jobs,
        "--pipe",
        "--recstart",
        "'>'",
        "--linebuffer",
        "prodigal",
        "-a",
        tmpprot_fpath,
        "-d",
        tmpnucl_fpath,
        "-q",
        "-p",
        "meta",
        "-m",
        "-o",
        os.devnull,
        "<",
        assembly,
        "2>",
        os.devnull,
    ]
    cmd = [str(arg) for arg in cmd]
    cmdline = subprocess.list2cmdline(cmd)
    logger.debug(cmdline)
    subprocess.run(cmdline, shell=True, check=True)
    search_path = os.path.join(tmpdir, "*.faa")
    aggregate_orfs(search_path, prots_out)
    search_path = os.path.join(tmpdir, "*.fna")
    aggregate_orfs(search_path, nucls_out)
    shutil.rmtree(tmpdir)


def run(
    assembly: str, nucls_out: str, prots_out: str, force: bool = False, cpus: int = 0
) -> Tuple[str, str]:
    """Calls ORFs from provided input assembly

    Parameters
    ----------
    assembly : str
        </path/to/assembly.fasta>
    nucls_out : str
        </path/to/nucls.out>
    prots_out : str
        </path/to/prots.out>
    force : bool
        overwrite outfpath if it already exists (the default is False).
    cpus : int
        num `cpus` to use. **Default (cpus=0) will run as many `cpus` as possible**

    Returns
    -------
    2-Tuple
        (`nucls_out`, `prots_out`)

    Raises
    -------
    FileExistsError
        `nucls_out` or `prots_out` already exists
    subprocess.CalledProcessError
        prodigal Failed
    ChildProcessError
        `nucls_out` or `prots_out` not written
    IOError
        `nucls_out` or `prots_out` incorrectly formatted
    """
    if not os.path.exists(assembly):
        raise FileNotFoundError(f"{assembly} does not exists!")

    if assembly.endswith(".gz"):
        assembly = gunzip(
            infpath=assembly, outfpath=assembly.rstrip(".gz"), delete_original=False
        )
    for fpath in [nucls_out, prots_out]:
        if os.path.exists(fpath) and not force:
            raise FileExistsError(f"{fpath} To overwrite use --force")
    if cpus == 1:
        annotate_sequential(assembly=assembly, prots_out=prots_out, nucls_out=nucls_out)
    else:
        annotate_parallel(
            assembly=assembly, prots_out=prots_out, nucls_out=nucls_out, cpus=cpus
        )
    for fp in [nucls_out, prots_out]:
        if not os.path.exists(fp) or not os.path.getsize(fp):
            raise ChildProcessError(f"{fp} not written")
        try:
            # Fasta file format check by simply reading orf annotations
            with open(fp) as fh:
                for _ in SimpleFastaParser(fh):
                    pass
        except (OSError, ValueError):
            raise OSError(f"{fp} file is not properly formatted")
    return nucls_out, prots_out


def contigs_from_headers(fpath: str) -> Mapping[str, str]:
    """Get ORF id to contig id translations using prodigal assigned ID from
    description.

    First determines if all of ID=3495691_2 from description is in header.
    "3495691_2" represents the 3,495,691st gene in the 2nd sequence.

    Example
    -------
    .. code-block:: python

        #: prodigal versions < 2.6 record
        >>>record.id
        'k119_1383959_3495691_2'

        >>>record.description
        'k119_1383959_3495691_2 # 688 # 1446 # 1 # ID=3495691_2;partial=01;start_type=ATG;rbs_motif=None;rbs_spacer=None'

        >>>record.description.split('#')[-1].split(';')[0].strip()
        'ID=3495691_2'

        >>>orf_id = '3495691_2'
        '3495691_2'

        >>>record.id.replace(f'_{orf_id}', '')
        'k119_1383959'

        #: prodigal versions >= 2.6 record
        >>>record.id
        'k119_1383959_2'
        >>>record.id.rsplit('_',1)[0]
        'k119_1383959'

    Parameters
    ----------
    fpath : str
        </path/to/prodigal/called/orfs.fasta>

    Returns
    -------
    dict
        contigs translated from prodigal ORF description.  {orf_id:contig_id, ...}

    """
    version = get_versions("prodigal")
    if version.count(".") >= 2:
        version = float(".".join(version.split(".")[:2]))
    else:
        version = float(version)
    translations = {}
    fh = gzip.open(fpath, "rt") if fpath.endswith(".gz") else open(fpath)
    for record in SeqIO.parse(fh, "fasta"):
        if version < 2.6:
            orf_id = (
                record.description.split("#")[-1]
                .split(";")[0]
                .strip()
                .replace("ID=", "")
            )
            contig_id = record.id.replace(f"_{orf_id}", "")
        else:
            contig_id = record.id.rsplit("_", 1)[0]
        translations.update({record.id: contig_id})
    fh.close()
    return translations


def orf_records_from_contigs(
    contigs: Union[List, Set], fpath: str
) -> List[SeqIO.SeqRecord]:
    """Retrieve list of *ORFs headers* from `contigs`. Prodigal annotated ORFs
    are required as the input `fpath`.

    Parameters
    ----------
    contigs: iterable
        iterable of contigs from which to retrieve ORFs
    fpath : str
        </path/to/prodigal/called/orfs.fasta>

    Returns
    -------
    list
        ORF SeqIO.SeqRecords from provided `contigs`. i.e. [SeqRecord, ...]

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    version = get_versions("prodigal")
    if version.count(".") >= 2:
        version = float(".".join(version.split(".")[:2]))
    else:
        version = float(version)

    records = []
    for record in SeqIO.parse(fpath, "fasta"):
        if version < 2.6:
            orf_id = (
                record.description.split("#")[-1]
                .split(";")[0]
                .strip()
                .replace("ID=", "")
            )
            contig_id = record.id.replace(f"_{orf_id}", "")
        else:
            contig_id = record.id.rsplit("_", 1)[0]
        if contig_id not in contigs:
            continue
        records.append(record)
    return records


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Calls ORFs with provided input assembly",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--assembly",
        help="Path to metagenome assembly",
        type=str,
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--output-nucls",
        help="Path to output nucleotide ORFs",
        type=str,
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--output-prots",
        help="Path to output amino-acid ORFs",
        type=str,
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--cpus",
        help="Number of processors to use. (If more than one this will parallelize prodigal using GNU parallel)",
        type=int,
        default=1,
        metavar="int",
    )
    parser.add_argument(
        "--force", help="Overwrite existing output ORF filepaths", action="store_true"
    )
    args = parser.parse_args()

    nucls_out, prots_out = run(
        assembly=args.assembly,
        nucls_out=args.output_nucls,
        prots_out=args.output_prots,
        cpus=args.cpus,
        force=args.force,
    )
    logger.info(f"written:\nnucls fpath: {nucls_out}\nprots fpath: {prots_out}")


if __name__ == "__main__":
    main()
