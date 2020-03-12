#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
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

Main script to run Autometa
"""


import logging
import os
import sys

import multiprocessing as mp

from autometa.config.user import AutometaUser
from autometa.common.utilities import timeit
from autometa.common.metagenome import Metagenome


logger = logging.getLogger('autometa')


def init_logger(fpath=None, level=None):
    """Initialize logger.

    By default will initialize streaming logger with DEBUG level messages.
    If `fpath` is provided, will write DEBUG level messages to `fpath` and
    set streaming messages to INFO.

    Parameters
    ----------
    fpath : str
        </path/to/file.log>
    level : int
        Overwrite default logging level behavior with provided `level`.
        This must be a constant from logging levels.
        See https://docs.python.org/3/library/logging.html#levels for details.
        i.e. logging.DEBUG, logging.INFO, etc. translates to 0,10, etc...

    Returns
    -------
    logging.Logger
        logging's Logger object to emit messages via methods:
        'warn','info','debug','error','exception','critical','fatal'

    Raises
    -------
    ValueError
        `level` must be int and one of 0, 10, 20, 30, 40, 50
    """
    levels = {
        logging.NOTSET,
        logging.DEBUG,
        logging.INFO,
        logging.WARNING,
        logging.ERROR,
        logging.CRITICAL}
    if level and type(level) is not int:
        raise ValueError(f'{level} must be an int! {type(level)}')
    if level and level not in levels:
        raise ValueError(f'{level} not in levels: {levels}!')
    formatter = logging.Formatter(
        fmt='[%(asctime)s %(levelname)s] %(name)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
    # Construct file/stream logging handlers
    streamhandler = logging.StreamHandler()
    streamhandler.setFormatter(formatter)
    if fpath:
        filehandler = logging.FileHandler(fpath)
        filehandler.setFormatter(formatter)
        logger.addHandler(filehandler)
        lvl = level if level else logging.INFO
    else:
        lvl = level if level else logging.DEBUG

    streamhandler.setLevel(lvl)
    logger.addHandler(streamhandler)
    logger.setLevel(logging.DEBUG)
    return logger

@timeit
def run(mgargs):
    """Run autometa.

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
    Subclassing an AutometaException class may be most appropriate use case here.
    """
    mg = Metagenome(
        assembly=mgargs.files.metagenome,
        outdir=mgargs.parameters.outdir,
        nucl_orfs_fpath=mgargs.files.nucleotide_orfs,
        prot_orfs_fpath=mgargs.files.amino_acid_orfs,
        taxonomy_fpath=mgargs.files.taxonomy,
        fwd_reads=mgargs.files.fwd_reads,
        rev_reads=mgargs.files.rev_reads,
        taxon_method=mgargs.parameters.taxon_method)
    try:
    # Original (raw) file should not be manipulated so return new object
        mg = mg.length_filter(
            out=mgargs.files.length_filtered,
            cutoff=mgargs.parameters.length_cutoff)
    except FileExistsError as err:
        logger.debug(f'{mgargs.files.length_filtered} already exists. Continuing..')
        mg = Metagenome(
            assembly=mgargs.files.length_filtered,
            outdir=mgargs.parameters.outdir,
            nucl_orfs_fpath=mgargs.files.nucleotide_orfs,
            prot_orfs_fpath=mgargs.files.amino_acid_orfs,
            taxonomy_fpath=mgargs.files.taxonomy,
            fwd_reads=mgargs.files.fwd_reads,
            rev_reads=mgargs.files.rev_reads,
            taxon_method=mgargs.parameters.taxon_method)
    # I.e. asynchronous execution here (work-queue tasks)
    mg.get_kmers(
        kmer_size=mgargs.parameters.kmer_size,
        normalized=mgargs.files.kmer_normalized,
        out=mgargs.files.kmer_counts,
        multiprocess=mgargs.parameters.kmer_multiprocess,
        nproc=mgargs.parameters.cpus,
        force=mgargs.parameters.force)

    coverages = mg.get_coverages(
        out=mgargs.files.coverages,
        from_spades=mgargs.parameters.cov_from_spades,
        sam=mgargs.files.sam,
        bam=mgargs.files.bam,
        lengths=mgargs.files.lengths,
        bed=mgargs.files.bed)
    # Filter by Kingdom
    kingdoms = mg.get_kingdoms(
        ncbi=mgargs.databases.ncbi,
        usepickle=mgargs.parameters.usepickle,
        blast=mgargs.files.blastp,
        hits=mgargs.files.blastp_hits,
        force=mgargs.parameters.force,
        cpus=mgargs.parameters.cpus)

    if not mgargs.parameters.kingdom in kingdoms:
        raise KeyError(f'{mgargs.parameters.kingdom} not recovered in dataset. Recovered: {", ".join(kingdoms.keys())}')
    mag = kingdoms.get(mgargs.parameters.kingdom)
    bins_df = mag.get_binning(
        method=mgargs.parameters.binning_method,
        kmers=mgargs.files.kmer_counts,
        embedded=mgargs.files.kmer_embedded,
        do_pca=mgargs.parameters.do_pca,
        pca_dims=mgargs.parameters.pca_dims,
        embedding_method=mgargs.parameters.embedding_method,
        coverage=coverages,
        domain=mgargs.parameters.kingdom,
        taxonomy=mgargs.files.taxonomy,
        reverse=mgargs.parameters.reversed,
    )
    binning_cols = ['cluster','completeness','purity']
    bins_df[binning_cols].to_csv(
        mgargs.files.binning,
        sep='\t',
        index=True,
        header=True)

def main(args):
    user = AutometaUser(dryrun=args.dryrun, nproc=args.cpus)
    for config in args.config:
        mgargs = user.prepare_run(config)
        run(mgargs)
        # cluster process -> mgargs.files.binning
        # TODO: Refine bins by connection mapping, taxon, or other methods
    # TODO: Construct pangenomes from multiple datasets
    # get_pangenomes()

if __name__ == '__main__':
    import argparse
    import time
    cpus = mp.cpu_count()
    parser = argparse.ArgumentParser('Main script to run Autometa')
    parser.add_argument('config',
        help='</path/to/metagenome.config>',
        nargs='*')
    parser.add_argument('--dryrun',
        help='whether to perform database updating/construction',
        action='store_true',
        default=False)
    parser.add_argument('--cpus',
        help=f'Num. cpus to use when updating/constructing databases (default: {cpus} cpus)',
        type=int,
        default=cpus)
    parser.add_argument('--debug',
        help=f'Stream debugging information to terminal',
        action='store_true',
        default=False)
    args = parser.parse_args()
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S",time.gmtime())
    level = logging.DEBUG if args.debug else None
    logger = init_logger(fpath=f'{timestamp}_autometa.log', level=level)
    try:
        main(args)
    except KeyboardInterrupt as err:
        logger.info('User cancelled run. Exiting...')
        sys.exit(1)
    except Exception as err:
        issue_request = '''
        An error was encountered!

        Please help us fix your problem!

        You may file an issue with us at https://github.com/KwanLab/Autometa/issues/new
        '''
        logger.exception(err)
        print(issue_request)
