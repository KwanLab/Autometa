#!/usr/bin/env python3
"""
Main script to run the Autometa pipeline
"""


import argparse
import os
import sys

import dependencies
from dataset import Dataset
from taxonomy import TaxonUtils


def main(args):
    dependencies.check_executables(args.verbose)
    config = dependencies.load_databases(args.config, args.verbose)
    for assembly in args.assemblies:
        for kingdom in ['bacteria','archaea']:
            dataset = Dataset(
                metagenome=assembly,
                dbconfig=config,
                marker_set='single_copy',
                kingdom=kingdom,
                usepickle=args.usepickle,
                verbose=args.verbose,
            )
            # Filter by length_cutoff
            dataset.filter_length(args.length_cutoff)
            # Call ORFs
            dataset.call_orfs()
            # Filter by Kingdom
            dataset.filter_taxonomy(
                rank_name=kingdom,
                by_rank='superkingdom',
                method='majority_vote')
            # # Get markers within filtered contigs to be used for bin assessment.
            dataset.get_markers()
            dataset.to_table()
            # # Bin using recursive dbscan or other method
            dataset.get_bins(method='recursive_dbscan')
            # # Refine bins by connection mapping, taxon, or other methods
            # mags.refine(by='connections')
            # mags.refine(by='taxa')
            # Construct pangenomes from multiple datasets
    # mags.get_pangenomes()

if __name__ == '__main__':
    parser = argparse.ArgumentParser('Script to retrieve taxonomy for Autometa pipeline')
    parser.add_argument('assemblies', nargs='+')
    parser.add_argument('--length-cutoff', default=3000, type=int)
    parser.add_argument(
        '--rank',
        default='superkingdom',
        choices=['superkingdom','phylum','class','order','family','genus','species'],
    )
    parser.add_argument('--method', default='majority_vote', choices=['majority_vote'])
    parser.add_argument('--outdir', default=os.curdir)
    parser.add_argument('--verbose', action='store_true', default=False)
    parser.add_argument('--force', action='store_true', default=False)
    parser.add_argument('--usepickle', action='store_true', default=False)
    parser.add_argument('--cpus',default=1, type=int)
    parser.add_argument(
        '--config',
        help='user defined databases config file',
        default=os.path.join(os.path.dirname(__file__), 'config', 'default.cfg'))
    args = parser.parse_args()
    main(args)
