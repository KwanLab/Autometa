#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Script to split metagenome assembly by kingdoms given the input votes. The lineages of the provided voted taxids will also be added and written to taxonomy.tsv
"""


import os
import logging
from pathlib import PurePath

import pandas as pd

from Bio import SeqIO
from typing import Union, List, Literal


from autometa.common.external import prodigal
from autometa.taxonomy import majority_vote
from autometa.taxonomy.lca import LCA
from autometa.taxonomy.ncbi import NCBI, NCBI_DIR
from autometa.taxonomy.gtdb import GTDB
from autometa.taxonomy.database import TaxonomyDatabase
from autometa.common.exceptions import TableFormatError

logger = logging.getLogger(__name__)


def assign(
    out: str,
    method: str = "majority_vote",
    assembly: str = None,
    prot_orfs: str = None,
    nucl_orfs: str = None,
    blast: str = None,
    lca_fpath: str = None,
    dbdir: str = NCBI_DIR,
    dbtype: Literal["ncbi", "gtdb"] = "ncbi",
    force: bool = False,
    verbose: bool = False,
    parallel: bool = True,
    cpus: int = 0,
) -> pd.DataFrame:
    """Assign taxonomy using `method` and write to `out`.

    Parameters
    ----------
    out : str
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
    lca_fpath : str, optional
        Path to output of LCA analysis, by default None
    dbdir : str, optional
        Path to NCBI databases directory, by default NCBI_DIR
    dbtype : str, optional
        Type of Taxonomy database to use, by default ncbi
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
    if os.path.exists(out) and os.path.getsize(out) and not force:
        logger.debug(f"FileExistsError: {out}. Use force to overwrite. skipping...")
        return pd.read_csv(out, sep="\t", index_col="contig")
    method = method.lower()
    if method != "majority_vote":
        raise NotImplementedError(method)
    outdir = os.path.dirname(os.path.realpath(out))
    lca_fpath = lca_fpath if lca_fpath else os.path.join(outdir, "lca.tsv")
    blast = blast if blast else os.path.join(outdir, "blastp.tsv")
    prot_orfs = prot_orfs if prot_orfs else os.path.join(outdir, "orfs.faa")
    nucl_orfs = nucl_orfs if nucl_orfs else os.path.join(outdir, "orfs.fna")

    taxa_db = NCBI(dbdir=dbdir) if dbtype == "ncbi" else GTDB(dbdir=dbdir)

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
            lca = LCA(taxonomy_db=taxa_db, verbose=verbose, cache=outdir)
        lca.blast2lca(
            blast=blast,
            out=lca_fpath,
            force=force,
        )

    def majority_vote_lca(out=out):
        if "lca" not in locals():
            lca = LCA(taxonomy_db=taxa_db, verbose=verbose, cache=outdir)
        ctg_lcas = lca.parse(lca_fpath=lca_fpath, orfs_fpath=prot_orfs)
        votes = majority_vote.rank_taxids(
            ctg_lcas=ctg_lcas, taxa_db=taxa_db, verbose=verbose
        )
        out = majority_vote.write_votes(results=votes, out=out)
        return pd.read_csv(out, sep="\t", index_col="contig")

    logger.info(f"Assigning taxonomy via {method}. This may take a while...")

    # Setup of taxonomy assignment sequence depending on file(s) provided
    calculation_sequence = {
        "lca_exists": [majority_vote_lca],
        "orfs_exists": [blast2lca, majority_vote_lca],
        "full": [call_orfs, blast2lca, majority_vote_lca],
    }
    # Now we need to determine which point to start the calculation...
    step = "full"
    for fp, argname in zip([lca_fpath, blast, prot_orfs], ["lca", "orfs", "orfs"]):
        if os.path.exists(fp) and os.path.getsize(fp):
            step = f"{argname}_exists"
            break

    if not assembly and step == "full":
        raise ValueError(f"assembly is required if no other files are specified!")

    logger.debug(f"starting taxonomy assignment sequence from {step}")
    for calculation in calculation_sequence[step]:
        logger.debug(f"running {calculation.__name__}")
        if calculation.__name__ == "majority_vote_lca":
            return calculation()
        calculation()


def add_ranks(df: pd.DataFrame, taxa_db: TaxonomyDatabase) -> pd.DataFrame:
    """Add canonical ranks to `df` and write to `out`

    Parameters
    ----------
    df : pd.DataFrame
        index="contig", column="taxid"
    taxa_db : TaxonomyDatabase
        NCBI or GTDB TaxonomyDatabase instance.

    Returns
    -------
    pd.DataFrame
        index="contig", columns=["taxid", *canonical_ranks]
    """
    dff = taxa_db.get_lineage_dataframe(df["taxid"].unique().tolist())
    return pd.merge(left=df, right=dff, how="left", left_on="taxid", right_index=True)


def get(
    filepath_or_dataframe: Union[str, pd.DataFrame],
    kingdom: str,
    taxa_db: TaxonomyDatabase,
) -> pd.DataFrame:
    """Retrieve specific `kingdom` voted taxa for `assembly` from `filepath`

    Parameters
    ----------
    filepath : str
        Path to tab-delimited taxonomy table. cols=['contig','taxid', *canonical_ranks]
    kingdom : str
        rank to retrieve from superkingdom column in taxonomy table.
    ncbi : str or autometa.taxonomy.NCBI instance, optional
        Path to NCBI database directory or NCBI instance, by default NCBI_DIR.
        This is necessary only if `filepath` does not already contain columns of canonical ranks.

    Returns
    -------
    pd.DataFrame
        DataFrame of contigs pertaining to retrieved `kingdom`.

    Raises
    ------
    FileNotFoundError
        Provided `filepath` does not exists or is empty.
    TableFormatError
        Provided `filepath` does not contain the 'superkingdom' column.
    KeyError
        `kingdom` is absent in provided taxonomy table.
    """
    if isinstance(filepath_or_dataframe, pd.DataFrame):
        df = filepath_or_dataframe
    elif isinstance(filepath_or_dataframe, str):
        if not os.path.exists(filepath_or_dataframe) or not os.path.getsize(
            filepath_or_dataframe
        ):
            raise FileNotFoundError(filepath_or_dataframe)
        df = pd.read_csv(filepath_or_dataframe, sep="\t", index_col="contig")
    elif isinstance(filepath_or_dataframe, PurePath):
        df = pd.read_csv(filepath_or_dataframe, sep="\t", index_col="contig")
    else:
        raise TypeError(f"{type(filepath_or_dataframe)}")
    if df.shape[1] <= 2:
        # Voting method will write out contig and its voted taxid (2 cols).
        # So here we add the canonical ranks using voted taxids.
        df = add_ranks(df=df, taxa_db=taxa_db)

    if "superkingdom" not in df.columns:
        raise TableFormatError(f"superkingdom is not in taxonomy columns {df.columns}")
    kingdom = kingdom.lower()
    df = df[df.superkingdom == kingdom]
    if df.empty:
        raise KeyError(f"{kingdom} not recovered in dataset.")
    return df


def write_ranks(
    taxonomy: pd.DataFrame,
    assembly: str,
    outdir: str,
    rank: str = "superkingdom",
    prefix: str = None,
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
    prefix : str, optional
        Prefix each of the paths written with `prefix` string.

    Returns
    -------
    list
        [rank_name_fpath, ...]

    Raises
    ------
    KeyError
        `rank` not in taxonomy columns
    """
    if rank not in taxonomy.columns:
        raise KeyError(f"{rank} not in taxonomy columns: {taxonomy.columns}")
    if not os.path.exists(assembly) or not os.path.getsize(assembly):
        raise FileNotFoundError(assembly)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    assembly_records = [record for record in SeqIO.parse(assembly, "fasta")]
    # Include unaligned records in unclassified fasta
    unaligned_contigs = set(
        record.id
        for record in SeqIO.parse(assembly, "fasta")
        if record.id not in taxonomy.index
    )
    fpaths = []
    for rank_name, dff in taxonomy.groupby(rank):
        # First determine the file path respective to the rank name
        rank_name = rank_name.replace(" ", "_")
        if prefix:
            rank_name_fname = ".".join([prefix, rank_name.lower(), "fna"])
        else:
            rank_name_fname = ".".join([rank_name.lower(), "fna"])
        rank_name_fpath = os.path.join(outdir, rank_name_fname)
        # Now retrieve and write records respective to rank
        # include unaligned contigs if rank is unclassified
        if rank_name == TaxonomyDatabase.UNCLASSIFIED:
            contig_set = set(dff.index).union(unaligned_contigs)
        else:
            contig_set = dff.index
        records = [record for record in assembly_records if record.id in contig_set]
        if not records:
            logger.warning(f"No records to write to {rank_name_fpath}")
        else:
            n_written = SeqIO.write(records, rank_name_fpath, "fasta")
            logger.debug(f"Wrote {n_written:,} records to {rank_name_fpath}")
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
    parser = argparse.ArgumentParser(
        description="Filter metagenome by taxonomy.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--votes",
        help="Input path to voted taxids table. should contain (at least) 'contig' and 'taxid' columns",
        type=str,
        metavar="filepath",
        required=True,
    )
    parser.add_argument(
        "--assembly",
        help="Path to metagenome assembly (nucleotide fasta).",
        metavar="filepath",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--output",
        help="Output directory to write specified canonical ranks fasta files and taxon-binning results table",
        type=str,
        metavar="dirpath",
        required=True,
    )
    parser.add_argument(
        "--prefix",
        help="prefix to use for each file written e.g. `prefix`.taxonomy.tsv. Note: Do not use a directory prefix.",
        type=str,
        metavar="str",
    )
    parser.add_argument(
        "--split-rank-and-write",
        help="If specified, will split contigs by provided canonical-rank column then write to `output` directory",
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
        "--dbdir",
        help="Path to taxonomy database directory.",
        metavar="dirpath",
        default=NCBI_DIR,
    )
    parser.add_argument(
        "--dbtype",
        help="Taxonomy database to use",
        choices=["ncbi", "gtdb"],
        default="ncbi",
    )
    args = parser.parse_args()

    if args.dbtype == "ncbi":
        taxa_db = NCBI(args.dbdir)
    elif args.dbtype == "gtdb":
        taxa_db = GTDB(args.dbdir)

    filename = f"{args.prefix}.taxonomy.tsv" if args.prefix else "taxonomy.tsv"
    out = os.path.join(args.output, filename)
    taxa_df = pd.read_csv(args.votes, sep="\t", index_col="contig")

    if not os.path.isdir(args.output):
        os.makedirs(args.output)

    if taxa_df.shape[1] <= 2:
        taxa_df = add_ranks(taxa_df, taxa_db=taxa_db)
        taxa_df.to_csv(out, sep="\t", index=True, header=True)
        logger.debug(
            f"Wrote {taxa_df.shape[0]:,} contigs canonical rank names to {out}"
        )
    if args.split_rank_and_write:
        if args.split_rank_and_write not in TaxonomyDatabase.CANONICAL_RANKS:
            raise ValueError(
                f"rank: {args.split_rank_and_write} not in {TaxonomyDatabase.CANONICAL_RANKS}"
            )
        written_ranks = write_ranks(
            taxonomy=taxa_df,
            assembly=args.assembly,
            outdir=args.output,
            rank=args.split_rank_and_write,
            prefix=args.prefix,
        )
        num_written_ranks = len(written_ranks)
        logger.info(f"Wrote {num_written_ranks} ranks to files.")


if __name__ == "__main__":
    main()
