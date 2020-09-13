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

Main script to run Autometa
"""


import logging
import os
import sys

import multiprocessing as mp
import pandas as pd

from autometa.common import coverage, kmers, markers, utilities
from autometa.common.metabin import MetaBin
from autometa.common.metagenome import Metagenome
from autometa.config.user import AutometaUser
from autometa.taxonomy import vote
from autometa.binning import recursive_dbscan
from autometa.common.exceptions import AutometaError


logger = logging.getLogger("autometa")


def init_logger(fpath=None, verbosity=0):
    """Initialize logger.

    By default will initialize streaming logger with INFO level messages.
    If `fpath` is provided, will write DEBUG level messages to `fpath` and
    set streaming messages to INFO.

    Parameters
    ----------
    fpath : str, optional
        </path/to/file.log>
    verbosity : int, optional
        Overwrite default logging level behavior with provided `verbosity`.
        This must be between 0-2 where 0 is the least information and 2 is the most information.
        See https://docs.python.org/3/library/logging.html#levels for details.

    Returns
    -------
    logging.Logger
        logging's Logger object to emit messages via methods:
        'warn','info','debug','error','exception','critical','fatal'

    Raises
    -------
    TypeError
        `verbosity` must be an int
    ValueError
        `verbosity` must be between of 0 and 2
    """
    log_levels = {
        0: logging.WARN,
        1: logging.INFO,
        2: logging.DEBUG,
    }

    if type(verbosity) is not int:
        raise TypeError(f"{verbosity} must be an int! {type(verbosity)}")
    if verbosity and verbosity not in log_levels:
        raise ValueError(f"{verbosity} not in log_levels: {log_levels}!")

    level = log_levels.get(verbosity)

    formatter = logging.Formatter(
        fmt="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    # Construct file/stream logging handlers
    streamhandler = logging.StreamHandler()
    streamhandler.setFormatter(formatter)
    if fpath:
        filehandler = logging.FileHandler(fpath)
        filehandler.setFormatter(formatter)
        logger.addHandler(filehandler)

    streamhandler.setLevel(level)
    logger.addHandler(streamhandler)
    logger.setLevel(logging.DEBUG)
    return logger


@utilities.timeit
def run_autometa(mgargs):
    """Run the autometa metagenome binning pipeline using the provided metagenome args.

    Pipeline
    --------

        #. Filter contigs in metagenome by length
        #. Determine metagenome coverages
        #. Count metagenome k-mers
        #. Call metagenome ORFs
        #. Assign contig taxonomy and filter by kingdom of interest
        #. Annotate gene markers
        #. Prepare primary table from annotations
        #. Bin contigs using primary table

    Parameters
    ----------
    mgargs : argparse.Namespace
        metagenome args

    Returns
    -------
    NoneType

    Raises
    -------
    TODO: Need to enumerate all exceptions raised from within binning pipeline.
    I.e. Demarkate new exception (not yet handled) vs. handled exception.
    Subclassing an AutometaError class may be most appropriate use case here.
    """

    # 1. Apply length filter (if desired).
    mg = Metagenome(
        assembly=mgargs.files.metagenome,
        outdir=mgargs.parameters.outdir,
        nucl_orfs_fpath=mgargs.files.nucleotide_orfs,
        prot_orfs_fpath=mgargs.files.amino_acid_orfs,
        fwd_reads=mgargs.files.fwd_reads,
        rev_reads=mgargs.files.rev_reads,
        se_reads=mgargs.files.se_reads,
    )
    # TODO asynchronous execution here (work-queue/Makeflow tasks) 1a.
    # Describe metagenome prior to length filter
    # COMBAK: Checkpoint length filtered
    try:
        # Original (raw) file should not be manipulated so we return a new object
        mg = mg.length_filter(
            out=mgargs.files.length_filtered, cutoff=mgargs.parameters.length_cutoff
        )
    except FileExistsError as err:
        logger.debug(f"{mgargs.files.length_filtered} already exists. Continuing..")
        mg = Metagenome(
            assembly=mgargs.files.length_filtered,
            outdir=mgargs.parameters.outdir,
            nucl_orfs_fpath=mgargs.files.nucleotide_orfs,
            prot_orfs_fpath=mgargs.files.amino_acid_orfs,
            fwd_reads=mgargs.files.fwd_reads,
            rev_reads=mgargs.files.rev_reads,
            se_reads=mgargs.files.se_reads,
        )
    # TODO asynchronous execution here (work-queue/Makeflow tasks) 1b.
    # Describe metagenome after length filter

    # TODO asynchronous execution here (work-queue/Makeflow tasks) 2a.
    counts = kmers.count(
        assembly=mg.assembly,
        size=mgargs.parameters.kmer_size,
        out=mgargs.files.kmer_counts,
        force=mgargs.parameters.force,
        verbose=mgargs.parameters.verbose,
        multiprocess=mgargs.parameters.kmer_multiprocess,
        cpus=mgargs.parameters.cpus,
    )
    kmers.normalize(
        df=counts,
        method=mgargs.parameters.kmer_transform,
        out=mgargs.files.kmer_normalized,
        force=mgargs.parameters.force,
    )
    # COMBAK: Checkpoint kmers
    # TODO asynchronous execution here (work-queue/Makeflow tasks) 2b.
    coverage.get(
        fasta=mg.assembly,
        out=mgargs.files.coverages,
        from_spades=mgargs.parameters.cov_from_spades,
        fwd_reads=mgargs.files.fwd_reads,
        rev_reads=mgargs.files.rev_reads,
        se_reads=mgargs.files.se_reads,
        sam=mgargs.files.sam,
        bam=mgargs.files.bam,
        lengths=mgargs.files.lengths,
        bed=mgargs.files.bed,
        cpus=mgargs.parameters.cpus,
    )
    # TODO asynchronous execution here (work-queue/Makeflow tasks) 2c.
    mg.call_orfs(
        force=mgargs.parameters.force,
        cpus=mgargs.parameters.cpus,
        parallel=mgargs.parameters.parallel,
    )

    # 3. Assign taxonomy (if desired).
    # COMBAK: Checkpoint coverages
    if mgargs.parameters.do_taxonomy:
        # First assign taxonomy
        vote.assign(
            method=mgargs.parameters.taxon_method,
            outfpath=mgargs.files.taxonomy,
            fasta=mg.assembly,
            prot_orfs=mgargs.files.amino_acid_orfs,
            nucl_orfs=mgargs.files.nucleotide_orfs,
            blast=mgargs.files.blastp,
            hits=mgargs.files.blastp_hits,
            lca_fpath=mgargs.files.lca,
            ncbi_dir=mgargs.databases.ncbi,
            tmpdir=mgargs.parameters.tmpdir,
            usepickle=mgargs.parameters.usepickle,
            force=mgargs.parameters.force,
            verbose=mgargs.parameters.verbose,
            parallel=mgargs.parameters.parallel,
            cpus=mgargs.parameters.cpus,
        )
        # Now filter by Kingdom
        metabin = vote.get(
            fpath=mgargs.files.taxonomy,
            assembly=mg.assembly,
            ncbi_dir=mgargs.databases.ncbi,
            kingdom=mgargs.parameters.kingdom,
            outdir=mgargs.parameters.outdir,
        )
        orfs_fpath = os.path.join(
            mgargs.parameters.outdir, f"{mgargs.parameters.kingdom}.orfs.faa"
        )
        metabin.write_orfs(fpath=orfs_fpath, orf_type="prot")
        contigs = set(metabin.contig_ids)
    else:
        orfs_fpath = mgargs.files.amino_acid_orfs
        contigs = {record.id for record in mg.seqrecords}

    # 4. Annotate marker genes for retained ORFs.
    # markers.get(...)
    # We use orfs_fpath here instead of our config filepath
    # because we may have filtered out ORFs by taxonomy.
    if mgargs.parameters.kingdom == "bacteria":
        scans = mgargs.files.bacteria_hmmscan
        markers_outfpath = mgargs.files.bacteria_markers
    else:
        scans = mgargs.files.archaea_hmmscan
        markers_outfpath = mgargs.files.archaea_markers

    markers_df = markers.get(
        kingdom=mgargs.parameters.kingdom,
        orfs=orfs_fpath,
        dbdir=mgargs.databases.markers,
        scans=scans,
        out=markers_outfpath,
        force=mgargs.parameters.force,
        seed=mgargs.parameters.seed,
    )

    # 5. Prepare contigs' annotations for binning
    # At this point, we should have kmer counts, coverages and possibly taxonomy info
    # Embed the kmer counts to retrieve our initial master dataframe
    master = kmers.embed(
        kmers=mgargs.files.kmer_normalized,
        out=mgargs.files.kmer_embedded,
        force=mgargs.parameters.force,
        embed_dimensions=mgargs.parameters.embed_dimensions,
        do_pca=mgargs.parameters.do_pca,
        pca_dimensions=mgargs.parameters.pca_dimensions,
        method=mgargs.parameters.embed_method,
        seed=mgargs.parameters.seed,
    )
    # Add coverage and possibly taxonomy annotations
    master = master[master.index.isin(contigs)]
    if mgargs.parameters.do_taxonomy:
        annotations = [mgargs.files.coverages, mgargs.files.taxonomy]
    else:
        annotations = [mgargs.files.coverages]
    for fpath in annotations:
        df = pd.read_csv(fpath, sep="\t", index_col="contig")
        master = pd.merge(master, df, how="left", left_index=True, right_index=True)
    master = master.convert_dtypes()

    # 6. Bin contigs
    bins_df = recursive_dbscan.binning(
        master=master,
        markers=markers_df,
        domain=mgargs.parameters.kingdom,
        completeness=mgargs.parameters.completeness,
        purity=mgargs.parameters.purity,
        taxonomy=mgargs.parameters.do_taxonomy,
        starting_rank=mgargs.parameters.starting_rank,
        method=mgargs.parameters.clustering_method,
        reverse_ranks=mgargs.parameters.reverse_ranks,
        verbose=mgargs.parameters.verbose,
    )
    binning_fpath = (
        mgargs.files.bacteria_binning
        if mgargs.parameters.kingdom == "bacteria"
        else mgargs.files.archaea_binning
    )
    binning_cols = ["cluster", "completeness", "purity"]
    bins_df[binning_cols].to_csv(binning_fpath, sep="\t", index=True, header=True)
    return bins_df


def main(args):
    """
    Main logic for running autometa pipeline.

    Warning: This should be called by `entrypoint` and not directly.

    Parameters
    ----------
    args : argparse.Namespace
        namespace containing config information for AutometaUser

    Returns
    -------
    NoneType
        Nothing if no errors are encountered.

    """

    logger = init_logger(fpath=args.log, verbosity=args.verbosity)
    # Configure AutometaUser
    # TODO: master from WorkQueue is AutometaUser
    user = AutometaUser(nproc=args.cpus)
    user.configure(dryrun=args.check_dependencies, update=args.update)

    # all_bins = {}
    for config in args.config:
        # TODO: Add directions to master from WorkQueue
        mgargs = user.prepare_binning_args(config)
        bins = run_autometa(mgargs)

        # TODO: refinement/processing/prep for pangenome algos
        # refined_bins = refine_binning(bins)
        # processed_bins = process_binning(refined_bins)
        # sample = mgargs.parameters.metagenome_num
        # all_bins.update({sample: processed_bins})
    # pangenomes = get_pangenomes(all_bins)


def entrypoint():
    """
    Main entrypoint for autometa pipeline.

    Note, a requirement of packaging and distribution is for entrypoints to not
    require any arguments. This is a wrapper to the main functionality of running
    autometa via the `main` function.

    Returns
    -------
    NoneType

    """
    import argparse

    # import time
    cpus = mp.cpu_count()
    parser = argparse.ArgumentParser(
        description="Main script to run the Autometa pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("config", help="Path to your metagenome.config file", nargs="*")
    parser.add_argument(
        "--cpus",
        help="Num. cpus to use when updating/constructing databases",
        type=int,
        default=cpus,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbosity",
        action="count",
        default=0,
        help="Verbosity (between 1-2 occurrences with more leading to more "
        "verbose logging). WARN=0, INFO=1, DEBUG=2",
    )
    parser.add_argument(
        "--log",
        help="Path to write a log file (e.g. </path/to/autometa.log>)",
        type=str,
    )
    parser.add_argument(
        "--check-dependencies",
        help="Check user executables and databases accessible to Autometa and exit.",
        action="store_true",
    )
    parser.add_argument(
        "--update",
        help="Update existing databases to most recent releases",
        action="store_true",
    )
    args = parser.parse_args()

    try:
        main(args)
    except KeyboardInterrupt:
        logger.info("User cancelled run. Exiting...")
        sys.exit(1)
    except Exception as err:
        logger.exception(err)
        if not issubclass(err.__class__, AutometaError):
            issue_request = """
            An error was encountered!

            Please help us fix your problem!

            You may file an issue with us at https://github.com/KwanLab/Autometa/issues/new/choose
            """
            logger.info(issue_request)
        sys.exit(1)


if __name__ == "__main__":
    entrypoint()
