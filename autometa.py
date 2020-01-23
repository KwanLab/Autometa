#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main script to run Autometa
"""


import logging
import os
import sys

from autometa.config.user import AutometaUser
from autometa.config import PROJECTS_DIR
from autometa.config import parse_config
from autometa.common.utilities import timeit
from autometa.common.metagenome import Metagenome


logger = logging.getLogger(__name__)

@timeit
def run(mgargs):
    """Run autometa.

    Parameters
    ----------
    mgargs : argparse.Namespace
        metagenome args

    Returns
    -------
    None
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

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
    binning_kingdoms = ['bacteria','archaea']
    kingdoms = mg.get_kingdoms(
        ncbi=mgargs.databases.ncbi,
        usepickle=mgargs.parameters.usepickle,
        blast=mgargs.files.blastp,
        hits=mgargs.files.blastp_hits,
        force=mgargs.parameters.force,
        cpus=mgargs.parameters.cpus)
    for kingdom,mag in kingdoms.items():
        # Get markers to be used for bin assessment.
        if kingdom not in binning_kingdoms:
            continue
        # markers = mag.markers(kingdom)
        bins_df = mag.get_binning(
            method=mgargs.parameters.binning_method,
            kmers=mgargs.files.kmer_counts,
            embedded=mgargs.files.kmer_embedded,
            do_pca=mgargs.parameters.do_pca,
            pca_dims=mgargs.parameters.pca_dims,
            embedding_method=mgargs.parameters.embedding_method,
            coverage=coverages,
            domain=kingdom,
            taxonomy=mgargs.files.taxonomy,
            reverse=mgargs.parameters.reversed,
        )
        binning_cols = ['cluster','completeness','purity']
        bins_df[binning_cols].to_csv(
            mgargs.files.binning,
            sep='\t',
            index=True,
            header=True)
        # TODO: Refine bins by connection mapping, taxon, or other methods
        # mag.refine(by='connections')
        # mag.refine(by='taxa')

def main(args):
    if not args.metagenomes_configs and not args.metagenomes and not args.resume:
        raise ValueError('Must provide metagenomes-configs or metagenomes')
    if args.config:
        user = AutometaUser(args.config, dryrun=args.dryrun)
    else:
        user = AutometaUser(dryrun=args.dryrun)
    # Configure environment and databases
    user.configure(nproc=args.cpus)
    # Workflow control...
    # TODO: WorkQueue handling. to process multiple metagenomes at once.
    if args.resume:
        mg_configs = user.get_mgargs(
            projects_dir=args.projects,
            project_num=args.project,
            metagenome_num=args.resume)
    elif args.metagenomes_configs:
        try:
            mg_configs = user.add_metagenomes(args.metagenomes_configs)
        except FileNotFoundError as err:
            project_configs = user.new_project(args)
            mg_configs = user.add_metagenomes(args.metagenomes_configs)
    else:
        project_configs = user.new_project(args)
        mg_configs =  project_configs.get('metagenomes')
    # Run autometa on workflow metagenome args...
    for metagenome,mgargs in mg_configs.items():
        run(mgargs)
    # user.bin_metagenome(metagenome_config)
    # TODO: Construct pangenomes from multiple datasets
    # get_pangenomes()

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)

    ###############################
    # AutometaUser Project(s) API #
    ###############################

    parser = argparse.ArgumentParser('Main script to run Autometa')
    parser.add_argument('--projects',
        help=f'</path/autometa/projects/dir> (Default is {PROJECTS_DIR}).',
        default=PROJECTS_DIR,
        required=False)
    parser.add_argument('--project',
        help='project number for which to resume autometa binning (required with `--resume` and --add-metagenome).',
        type=int)
    parser.add_argument('--resume',
        help='metagenome number for which to resume autometa binning (`--project` num is required).',
        type=int,
        default=0)
    parser.add_argument('--metagenomes-configs',
        help='</path/to/metagenome.config>',
        nargs='*')
    parser.add_argument('--dryrun',
        help='whether to perform database updating/construction',
        action='store_true',
        default=False)

    #######################
    # Autometa Parameters #
    #######################

    parser.add_argument('metagenomes', nargs='*')
    parser.add_argument('--length-cutoff', default=3000, type=int)
    parser.add_argument('--cov-from-spades',
        help='retrieve coverage from spades headers. (Only may be used when SPAdes assemblies are provided)',
        action='store_true',
        default=False)
    parser.add_argument(
        '--kmer-size',
        help='size of k-mer to calculate frequencies.',
        default=5, type=int)
    parser.add_argument(
        '--kmer-multiprocess',
        help='use multiprocessing to count k-mers.',
        action='store_true', default=False)
    parser.add_argument(
        '--kmer-normalize',
        help='Perform CLR transform on k-mer frequencies.',
        action='store_true', default=False)
    parser.add_argument('--do-pca',
        help='Perform PCA prior to running embedding method', default=False, action='store_true')
    parser.add_argument(
        '--pca-dims',
        help='Number of dimesions to reduce k-mer frequencies using PCA',
        default=50, type=int)
    parser.add_argument(
        '--embedding-method',
        help='Embedding method for dimension reduction of contig k-mer frequencies',
        default='UMAP',
        choices=['TSNE','UMAP'])
    parser.add_argument('--taxon-method', default='majority_vote', choices=['majority_vote'])
    parser.add_argument('--reversed', help='Reverse order at which taxonomic ranks are clustered', default=True, action='store_false')
    parser.add_argument('--binning-method',
        default='recursive_dbscan',
        choices=['recursive_dbscan'])
    parser.add_argument('--completeness', type=float, default=20.)
    parser.add_argument('--purity', type=float, default=90.)
    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--force', action='store_true', default=False)
    parser.add_argument('--usepickle', action='store_true', default=False)
    parser.add_argument('--parallel', help="Use GNU parallel",
        action='store_true', default=False)
    parser.add_argument('--cpus',default=1, type=int)
    parser.add_argument('--config',help='user defined config file')
    args = parser.parse_args()
    main(args)
