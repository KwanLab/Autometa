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

taxonomy voting script
"""


import os
import logging
import tempfile
import shutil

import pandas as pd

from Bio import SeqIO
from typing import Union, List


from autometa.common.external import prodigal
from autometa.common.metabin import MetaBin
from autometa.taxonomy import majority_vote
from autometa.taxonomy.lca import LCA
from autometa.taxonomy.ncbi import NCBI_DIR
from autometa.taxonomy.ncbi import NCBI
from autometa.common.exceptions import TableFormatError

logger = logging.getLogger(__name__)


def assign(
    outfpath: str,
    method: str = "majority_vote",
    assembly: str = None,
    prot_orfs: str = None,
    nucl_orfs: str = None,
    blast: str = None,
    hits: str = None,
    lca_fpath: str = None,
    ncbi_dir: str = NCBI_DIR,
    tmpdir: str = None,
    usepickle: bool = True,
    force: bool = False,
    verbose: bool = False,
    parallel: bool = True,
    cpus: int = 0,
) -> pd.DataFrame:
    """Assign taxonomy using `method` and write to `outfpath`.

    Parameters
    ----------
    outfpath : str
        Path to write taxonomy table of votes
    method : str, optional
        Method to assign contig taxonomy, by default "majority_vote".
        choices include "majority_vote", ...
    assembly : str, optional
        Path to assembly fasta file (nucleotide), by default None
    prot_orfs : str, optional
        Path to amino-acid ORFs called from `assembly`, by default None
    nucl_orfs : str, optional
        Path to nucleotide ORFs called from `assembly`, by default None
    blast : str, optional
        Path to blastp table, by default None
    hits : str, optional
        Path to pickled hits parsed from blastp table, by default None
    lca_fpath : str, optional
        Path to output of LCA analysis, by default None
    ncbi_dir : str, optional
        Path to NCBI databases directory, by default NCBI_DIR
    tmpdir : str, optional
        Path to use as temporary directory, by default None
    usepickle : bool, optional
        Pickle LCA analysis data structures for later LCA queries, by default True
    force : bool, optional
        Overwrite existing annotations, by default False
    verbose : bool, optional
        Increase verbosity, by default False
    parallel : bool, optional
        Whether to perform annotations using multiprocessing and GNU parallel, by default True
    cpus : int, optional
        Number of cpus to use if `parallel` is True, by default will try to use all available.

    Returns
    -------
    pd.DataFrame
        index="contig", columns=["taxid"]

    Raises
    ------
    NotImplementedError
        Provided `method` has not yet been implemented.
    ValueError
        Assembly file is required if no other annotations are provided.
    """
    if os.path.exists(outfpath) and os.path.getsize(outfpath) and not force:
        logger.debug(
            f"FileExistsError: {outfpath}. Use force to overwrite. skipping..."
        )
        return pd.read_csv(outfpath, sep="\t", index_col="contig")
    method = method.lower()
    if method != "majority_vote":
        raise NotImplementedError(method)
    outdir = os.path.dirname(os.path.realpath(outfpath))
    tmpdir = (
        tmpdir
        if tmpdir
        else tempfile.mkdtemp(suffix=None, prefix="taxon-assignment", dir=outdir)
    )
    lca_fpath = lca_fpath if lca_fpath else os.path.join(tmpdir, "lca.tsv")
    hits = hits if hits else os.path.join(tmpdir, "hits.pkl.gz")
    blast = blast if blast else os.path.join(tmpdir, "blastp.tsv")

    def call_orfs():
        prodigal.run(
            assembly=assembly,
            nucls_out=nucl_orfs,
            prots_out=prot_orfs,
            force=force,
            cpus=cpus,
            parallel=parallel,
        )

    def blast2lca():
        if "lca" not in locals():
            lca = LCA(
                dbdir=ncbi_dir,
                outdir=outdir,
                usepickle=usepickle,
                verbose=verbose,
                cpus=cpus,
            )
        lca.blast2lca(
            fasta=prot_orfs,
            outfpath=lca_fpath,
            blast=blast,
            hits_fpath=hits,
            force=force,
        )

    def majority_vote_lca():
        if "lca" not in locals():
            lca = LCA(
                dbdir=ncbi_dir,
                outdir=outdir,
                usepickle=usepickle,
                verbose=verbose,
                cpus=cpus,
            )
        ctg_lcas = lca.parse(lca_fpath=lca_fpath, orfs_fpath=prot_orfs)
        votes = majority_vote.rank_taxids(ctg_lcas=ctg_lcas, ncbi=lca)
        outfpath = majority_vote.write_votes(results=votes, outfpath=outfpath)
        return pd.read_csv(outfpath, sep="\t", index_col="contig")

    logger.info(f"Assigning taxonomy via {method}. This may take a while...")

    # Setup of taxonomy assignment sequence depending on file(s) provided
    calculation_sequence = {
        "lca_exists": [majority_vote_lca],
        "orfs_exists": [blast2lca, majority_vote_lca],
        "full": [call_orfs, blast2lca, majority_vote_lca],
    }
    # Now we need to determine which point to start the calculation...
    step = "full"
    for fp, argname in zip(
        [lca_fpath, hits, blast, prot_orfs], ["lca", "lca", "lca", "orfs"],
    ):
        if os.path.exists(fp):
            step = f"{argname}_exists"
            break

    if assembly and step == "full":
        raise ValueError(f"assembly is required if no other files are specified!")

    logger.debug(f"starting taxonomy assignment sequence from {step}")
    try:
        for calculation in calculation_sequence[step]:
            logger.debug(f"running {calculation.__name__}")
            if calculation.__name__ == "majority_vote_lca":
                return calculation()
            calculation()
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def add_ranks(df: pd.DataFrame, outfpath: str, ncbi: Union[NCBI, str]) -> pd.DataFrame:
    """Add canonical ranks to `df` and write to `outfpath`

    Parameters
    ----------
    df : pd.DataFrame
        index="contig", column="taxid"
    outfpath : str
        Path to write taxonomy table with canonical ranks added.
    ncbi : str, NCBI
        Path to NCBI databases directory, or autometa NCBI instance.

    Returns
    -------
    pd.DataFrame
        index="contig", columns=["taxid", *canonical_ranks]
    """
    ncbi = ncbi if isinstance(ncbi, NCBI) else NCBI(ncbi)
    dff = ncbi.get_lineage_dataframe(df["taxid"].unique().tolist())
    df = pd.merge(left=df, right=dff, how="left", left_on="taxid", right_index=True,)
    # This overwrites the existing table with the canonical ranks added.
    df.to_csv(outfpath, sep="\t", index=True, header=True)
    # COMBAK: Add checkpointing checksum check here
    logger.debug(f"Added canonical rank names to {outfpath}")
    return df


def get(
    fpath: str,
    assembly: str,
    kingdom: str,
    ncbi_dir: str = NCBI_DIR,
    outdir: str = None,
) -> MetaBin:
    """Retrieve specific `kingdom` voted taxa for `assembly` from `fpath`

    Parameters
    ----------
    fpath : str
        Path to tab-delimited taxonomy table. cols=['contig','taxid', *canonical_ranks]
    assembly : str
        Path to assembly fasta file
    kingdom : str
        rank to retrieve from superkingdom column in taxonomy table.
    ncbi_dir : str, optional
        Path to NCBI database directory, by default NCBI_DIR.
        This is necessary only if `fpath` does not already contain columns of canonical ranks.
    outdir : str, optional
        Path to output directory, by default None

    Returns
    -------
    autometa.common.MetaBin
        MetaBin of contigs pertaining to retrieved `kingdom`.

    Raises
    ------
    FileNotFoundError
        Provided `fpath` does not exists or is empty.
    TableFormatError
        Provided `fpath` does not contain the 'superkingdom' column.
    KeyError
        `kingdom` is absent in provided taxonomy table.
    """
    if not os.path.exists(fpath) or not os.path.getsize(fpath):
        raise FileNotFoundError(fpath)
    df = pd.read_csv(fpath, sep="\t", index_col="contig")
    if df[1] <= 2:
        # Voting method will write out contig and its voted taxid (2 cols).
        # So here we add the canonical ranks using voted taxids.
        df = add_ranks(df=df, outfpath=fpath, ncbi_dir=ncbi_dir)

    if "superkingdom" not in df.columns:
        raise TableFormatError(f"superkingdom is not in taxonomy columns {df.columns}")
    kingdom = kingdom.lower()
    df = df[df.superkingdom == kingdom]
    if df.empty:
        raise KeyError(f"{kingdom} not recovered in dataset.")
    return MetaBin(assembly=assembly, contig_ids=df.index.tolist(), outdir=outdir)


def write_ranks(
    taxonomy: pd.DataFrame, assembly: str, outdir: str, rank: str = "superkingdom"
) -> List[str]:
    """Write fastas split by `rank`

    Parameters
    ----------
    taxonomy : pd.DataFrame
        dataframe containing canonical ranks of contigs assigned from :func:autometa.taxonomy.vote.assign(...)
    assembly : str
        Path to assembly fasta file
    outdir : str
        Path to output directory to write fasta files
    rank : str, optional
        canonical rank column in taxonomy table to split by, by default "superkingdom"

    Returns
    -------
    list
        [rank_name_fpath, ...]

    Raises
    ------
    ValueError
        `rank` not in canonical ranks
    """
    if rank not in NCBI.CANONICAL_RANKS:
        raise ValueError(f"rank: {rank} not in {NCBI.CANONICAL_RANKS}")

    assembly_records = [r for r in SeqIO.parse(assembly, "fasta")]
    fpaths = []
    for rank_name, dff in taxonomy.groupby(rank):
        # First determine the file path respective to the rank name
        rank_name = rank_name.replace(" ", "_")
        rank_name_fname = ".".join([rank_name.title(), "fna"])
        rank_name_fpath = os.path.join(outdir, rank_name_fname)
        # Now retrieve and write records respective to rank
        records = [record for record in assembly_records if record.id in dff.index]
        if not records:
            logger.warning(f"No records to write to {rank_name_fpath}")
        else:
            n_written = SeqIO.write(records, rank_name_fpath, "fasta")
            logger.debug(f"Wrote {n_written} records to {rank_name_fpath}")
            fpaths.append(rank_name_fpath)
    # Finally we will return all of the file paths that were written
    return fpaths


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    # Note: If you do not have any defaults corresponding to your parameters,
    # you may remove the formatter class: ArgumentDefaultsHelpFormatter
    # to reduce help text verbosity.
    parser = argparse.ArgumentParser(
        description="""
    This script handles filtering by taxonomy and can calculate various metagenome statistics.
    """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "taxonomy", help="Output path to write taxonomy table", type=str
    )
    parser.add_argument(
        "outdir", help="Output directory to store annotations.", type=str
    )
    parser.add_argument(
        "--assembly", help="Path to metagenome assembly (nucleotide fasta).", type=str,
    )
    parser.add_argument(
        "--nucls",
        help="Path to nucleotide ORFs corresponding to `assembly.` "
        "(Will write to path if ORFs do not exist).",
        type=str,
    )
    parser.add_argument(
        "--prots",
        help="Path to amino acid ORFs corresponding to `assembly.` "
        "(Will write to path if ORFs do not exist).",
        type=str,
    )
    parser.add_argument(
        "--split-rank-and-write",
        help="If specified, will split contigs by provided canonical-rank column then write to `outdir`",
        choices=[
            "superkingdom",
            "phylum",
            "class",
            "order",
            "family",
            "genus",
            "species",
        ],
    )
    parser.add_argument(
        "--taxon-method", default="majority_vote", choices=["majority_vote"]
    )
    parser.add_argument(
        "--ncbi", help="Path to NCBI databases directory.", default=NCBI_DIR
    )
    parser.add_argument(
        "--pickle",
        help="Whether to serialize taxonomy-specific files",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--blast",
        help="Path to diamond blast results.tsv. (outfmt=6) "
        "(Will write to path if it does not exist).",
        default=None,
    )
    parser.add_argument(
        "--hits",
        help="Path to diamond blast hits.pkl.gz. "
        "(Will write to path if it does not exist and `--pickle` is specified).",
        default=None,
    )
    # Eventually will need to create a subparser for the taxon assignment methods
    # to include help information and required parameters according to method.
    parser.add_argument(
        "--cpus",
        help="Number of cpus to use when performing "
        "annotations (Default will use all possible if `--parallel` is supplied).",
        default=0,
        type=int,
    )
    parser.add_argument(
        "--parallel",
        help="Use GNU parallel when performing annotations.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--force", help="Overwrite existing files", action="store_true", default=False
    )
    parser.add_argument(
        "--verbose",
        help="Log more information to terminal.",
        action="store_true",
        default=False,
    )

    args = parser.parse_args()

    assign(
        method=args.method,
        outfpath=args.taxonomy,
        fasta=args.assembly,
        prot_orfs=args.prot_orfs,
        nucl_orfs=args.nucl_orfs,
        blast=args.blast,
        hits=args.hits,
        lca_fpath=args.lca_fpath,
        ncbi_dir=args.ncbi,
        tmpdir=args.tmpdir,
        usepickle=args.usepickle,
        force=args.force,
        verbose=args.verbose,
        parallel=args.parallel,
        cpus=args.cpus,
    )
    get(
        fpath=args.taxonomy,
        assembly=args.assembly,
        kingdom=args.kingdom,
        ncbi_dir=args.ncbi,
        outdir=args.outdir,
    )
    written_ranks = write_ranks(
        taxonomy=args.taxonomy,
        assembly=args.assembly,
        outdir=args.outdir,
        rank=args.split_rank_and_write,
    )
    num_written_ranks = len(written_ranks)
    logger.info(f"Wrote {num_written_ranks} ranks to files.")


if __name__ == "__main__":
    main()
