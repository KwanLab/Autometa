#!/usr/bin/env python
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

This script contains the modified majority vote algorithm used in Autometa version 1.0
"""


import logging
import os
import sys
from typing import Dict, Union

from tqdm import tqdm

from autometa.taxonomy.lca import LCA
from autometa.taxonomy.ncbi import NCBI, NCBI_DIR

logger = logging.getLogger(__name__)


def is_consistent_with_other_orfs(
    taxid: int, rank: str, rank_counts: Dict[str, Dict], ncbi: NCBI
) -> bool:
    """Determines whether the majority of proteins in a contig, with rank equal
    to or above the given rank, are common ancestors of the taxid.

    If the majority are, this function returns True, otherwise it returns False.

    Parameters
    ----------
    taxid : int
        `taxid` to search against other taxids at `rank` in `rank_counts`.
    rank : str
        Canonical rank to search in `rank_counts`.
        Choices: species, genus, family, order, class, phylum, superkingdom.
    rank_counts : dict
        LCA canonical rank counts retrieved from ORFs respective to a contig.
        e.g. {canonical_rank: {taxid: num_hits, ...}, ...}
    ncbi : NCBI instance
        Instance or subclass of NCBI from autometa.taxonomy.ncbi.

    Returns
    -------
    boolean
        If the majority of ORFs in a contig are equal or above given rank then
        return True, otherwise return False.

    """
    rank_index = NCBI.CANONICAL_RANKS.index(rank)
    ranks_to_consider = NCBI.CANONICAL_RANKS[rank_index:]
    # Now we total up the consistent and inconsistent ORFs
    consistent = 0
    inconsistent = 0
    for rank_name in ranks_to_consider:
        if rank_name not in rank_counts:
            continue
        for rank_taxid, count in rank_counts[rank_name].items():
            if ncbi.is_common_ancestor(rank_taxid, taxid):
                consistent += count
            else:
                inconsistent += count
    if consistent > inconsistent:
        # COMBAK: See issue-#48: This could also return the ratio of consistent
        # to inconsistent to give the user an idea of the consistency of the
        # taxon assignments.
        return True
    else:
        return False


def lowest_majority(rank_counts: Dict[str, Dict], ncbi: NCBI) -> int:
    """Determine the lowest majority given `rank_counts` by first attempting to
    get a taxid that leads in counts with the highest specificity in terms of
    canonical rank.

    Parameters
    ----------
    rank_counts : dict
        {canonical_rank:{taxid:num_hits, ...}, rank2: {...}, ...}
    ncbi : NCBI instance
        NCBI object from autometa.taxonomy.ncbi

    Returns
    -------
    int
        Taxid above the lowest majority threshold.

    """
    taxid_totals = {}
    for rank in NCBI.CANONICAL_RANKS:
        if rank not in rank_counts:
            continue
        rank_index = NCBI.CANONICAL_RANKS.index(rank)
        ranks_to_consider = NCBI.CANONICAL_RANKS[rank_index:]
        for taxid in rank_counts[rank]:
            # Make a dictionary to total the number of canonical ranks hit
            # while traversing the path so that we can add 'unclassified' to
            # any that don't exist. Later we need to make sure that
            # 'unclassified' doesn't ever win
            ranks_in_path = {
                rank_to_consider: 0 for rank_to_consider in ranks_to_consider
            }
            # We need to add to taxid_totals for each taxid in the tax_path
            current_taxid = taxid
            current_rank = rank
            while current_taxid != 1:
                if current_rank not in set(NCBI.CANONICAL_RANKS):
                    current_taxid = ncbi.parent(current_taxid)
                    current_rank = ncbi.rank(current_taxid)
                    continue
                ranks_in_path[current_rank] += 1
                if current_rank not in taxid_totals:
                    taxid_totals.update({current_rank: {current_taxid: 1}})
                    current_taxid = ncbi.parent(current_taxid)
                    current_rank = ncbi.rank(current_taxid)
                    continue
                if current_taxid in taxid_totals[current_rank]:
                    taxid_totals[current_rank][current_taxid] += 1
                else:
                    taxid_totals[current_rank][current_taxid] = 1
                current_taxid = ncbi.parent(current_taxid)
                current_rank = ncbi.rank(current_taxid)
            # Now go through ranks_in_path. Where total = 0, add 'unclassified'
            for rank_to_consider in ranks_to_consider:
                if ranks_in_path[rank_to_consider] == 0:
                    if rank_to_consider not in taxid_totals:
                        taxid_totals[rank_to_consider] = {"unclassified": 1}
                    elif "unclassified" in taxid_totals[rank_to_consider]:
                        taxid_totals[rank_to_consider]["unclassified"] += 1
                    else:
                        taxid_totals[rank_to_consider]["unclassified"] = 1
    # If there are any gaps in the taxonomy paths for any of the proteins in the contig,
    # we need to add 'unclassified' to the relevant canonical taxonomic rank.
    # However, we must never allow 'unclassified' to win! (That just won't really tell us anything)
    # Now we need to determine which is the first level to have a majority
    for rank in NCBI.CANONICAL_RANKS:
        total_votes = 0
        taxid_leader = None
        taxid_leader_votes = 0
        if not rank in taxid_totals:
            continue
        for taxid in taxid_totals[rank]:
            taxid_votes = taxid_totals[rank][taxid]
            total_votes += taxid_votes
            if taxid_votes > taxid_leader_votes:
                taxid_leader = taxid
                taxid_leader_votes = taxid_votes
        majority_threshold = float(total_votes) / 2
        if taxid_leader_votes > majority_threshold and taxid_leader != "unclassified":
            return taxid_leader
    # Just in case
    return 1


def rank_taxids(
    ctg_lcas: dict, ncbi: Union[NCBI, LCA], verbose: bool = False
) -> Dict[str, int]:
    """Votes for taxids based on modified majority vote system where if a
    majority does not exist, the lowest majority is voted.

    Parameters
    ----------
    ctg_lcas : dict
        {ctg1:{canonical_rank:{taxid:num_hits,...},...}, ctg2:{...},...}
    ncbi : ncbi.NCBI or lca.LCA object
        instance of NCBI subclass or NCBI containing NCBI methods.
    verbose : bool
        Description of parameter `verbose` (the default is False).

    Returns
    -------
    dict
        {contig:voted_taxid, contig:voted_taxid, ...}

    """
    logging.debug("Ranking taxids")
    n_contigs = len(ctg_lcas) if verbose else None
    disable = False if verbose else True
    desc = "Ranking taxids" if verbose else None
    top_taxids = {}
    for contig in tqdm(
        ctg_lcas, disable=disable, total=n_contigs, desc=desc, leave=False
    ):
        acceptedTaxid = None
        for rank in NCBI.CANONICAL_RANKS:
            if acceptedTaxid is not None:
                break
            # Order in descending order of votes
            if rank in ctg_lcas[contig]:
                ordered_taxids = sorted(
                    ctg_lcas[contig][rank],
                    key=lambda tid: ctg_lcas[contig][rank][tid],
                    reverse=True,
                )
                for taxid in ordered_taxids:
                    if is_consistent_with_other_orfs(
                        taxid, rank, ctg_lcas[contig], ncbi
                    ):
                        acceptedTaxid = taxid
                        break
        # If acceptedTaxid is still None at this point, there was some kind of
        # draw, so we need to find the lowest taxonomic level where there is a
        # majority
        if acceptedTaxid is None:
            acceptedTaxid = lowest_majority(ctg_lcas[contig], ncbi)
        top_taxids[contig] = acceptedTaxid
    return top_taxids


def write_votes(results: Dict[str, int], out: str) -> str:
    """Writes voting `results` to provided `outfpath`.

    Parameters
    ----------
    results : dict
        {contig:voted_taxid, contig:voted_taxid, ...}
    out : str
        </path/to/results.tsv>.

    Returns
    -------
    str
        </path/to/results.tsv>

    Raises
    -------
    FileExistsError
        Voting results file already exists

    """
    if os.path.exists(out) and os.path.getsize(out):
        raise FileExistsError(out)
    lines = "contig\ttaxid\n"
    fh = open(out, "w")
    for contig, taxid in results.items():
        if sys.getsizeof(lines) >= 60000:
            fh.write(lines)
            lines = ""
        lines += f"{contig}\t{taxid}\n"
    fh.write(lines)
    fh.close()
    return out


def majority_vote(
    orfs: str,
    out: str,
    ncbi_dir: str,
    lca_out: str = None,
    verbose: bool = False,
    blast: str = None,
    force: bool = False,
) -> str:
    """Wrapper for modified majority voting algorithm from Autometa 1.0

    Parameters
    ----------
    orfs : str
        Path to ORFs fasta containing amino-acid sequences to be annotated.
    out : str
        Path to write assigned taxids.
    ncbi_dir : str
        Path to NCBI databases directory.
    lca_out : str, optional
        Path to write lowest common ancestor assignments table.
    verbose : bool, optional
        Increase verbosity of logging stream
    blast : str, optional
        Path to blast table (Note: Must be outfmt 6).
    force : bool, optional
        Whether to overwrite existing LCA results.

    Returns
    -------
    str
        Path to assigned taxids table.

    """
    lca = LCA(dbdir=ncbi_dir, verbose=verbose)
    if not lca_out:
        filename, __ = os.path.splitext(os.path.basename(orfs))
        outdir = os.path.dirname(os.path.realpath(out))
        lca_out = os.path.join(outdir, ".".join([filename, "lca.tsv"]))
    lca_fpath = lca.blast2lca(
        orfs=orfs,
        out=lca_out,
        blast=blast,
        force=force,
    )
    # retrieve lca taxids for each contig
    classifications = lca.parse(lca_fpath=lca_fpath, orfs_fpath=orfs)
    # Vote for majority lca taxid from contig lca taxids
    # We can pass in lca as the ncbi object here, because it is a subclass of NCBI.
    voted_taxids = rank_taxids(ctg_lcas=classifications, ncbi=lca, verbose=verbose)
    return write_votes(results=voted_taxids, out=out)


def main():
    import argparse

    basedir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
    parser = argparse.ArgumentParser(
        description="Script to assign taxonomy via a modified majority voting"
        " algorithm.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--orfs",
        help="Path to ORFs fasta containing amino-acid sequences to be annotated.",
        required=True,
    )
    parser.add_argument(
        "--output", help="Path to write voted taxid results table.", required=True
    )
    parser.add_argument(
        "--dbdir", help="Path to NCBI databases directory.", default=NCBI_DIR
    )
    parser.add_argument("--lca", help="Path to LCA results table.")
    parser.add_argument(
        "--blast",
        help="Path to BLASTP results table. "
        "Note: Results must be formatted using outfmt=6.",
    )
    parser.add_argument(
        "--verbose",
        help="Add verbosity to logging stream.",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    majority_vote(
        orfs=args.orfs,
        out=args.output,
        ncbi_dir=args.dbdir,
        blast=args.blast,
        lca_out=args.lca,
        verbose=args.verbose,
    )


if __name__ == "__main__":
    main()
