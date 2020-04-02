#!/usr/bin/env python
"""
COPYRIGHT
Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan

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

Processes metagenome assembled genomes from autometa binning results
"""


import argparse
import logging
import os
import subprocess

import pandas as pd

from Bio import SeqIO

logger = logging.getLogger(__name__)


def run_command(command_string, stdout_path = None):
    # Function that checks if a command ran properly. If it didn't, then log an error message then quit
    logger.debug('cluster_process.py, run_command: ' + command_string)
    if stdout_path:
        f = open(stdout_path, 'w')
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        logger.error('cluster_process.py: Error, the command:')
        logger.error(command_string)
        logger.error('failed, with exit code ' + str(exit_code))
        exit(1)

def assess_assembly(seq_record_list):
    assembly_size = sum(len(seq) for seq in seq_record_list)
    number_of_sequences = len(seq_record_list)
    sorted_seqs = sorted(seq_record_list, key=len)
    largest_sequence_length = len(sorted_seqs[-1])
    sequence_total = 0
    n50 = None
    for i,seq_object in enumerate(sorted_seqs):
        sequence_total += len(seq_object)
        if sequence_total > (float(assembly_size)/2):
            n50 = len(seq_object)
            break
    return {'size': assembly_size,
        'number_sequences':number_of_sequences,
        'largest_sequence':largest_sequence_length,
        'n50':n50}

def main(args):
    # Check paths exist
    if not os.path.isfile(args.bin_table):
        logger.error(f'Error! Could not find a bin table at the following path: {args.bin_table}')
        exit(1)

    if not os.path.isfile(args.fasta):
        logger.error('Error! Cannot find a fasta file at the following path: ' + args.fasta)
        exit(1)

    # Make output directory if it isn't already there
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    master_table = pd.read_csv(args.bin_table, sep='\t')

    # Format check for the table
    columns_to_check = [args.column, 'contig', 'length', 'cov', 'single_copy_PFAMs']
    if args.do_taxonomy:
        columns_to_check.append('taxid')

    for column in columns_to_check:
        if column not in master_table.columns:
            logger.error('Error! Could not find a column called ' + column + ' in table ' + args.bin_table)
            exit(1)

    contig_info = dict() # Will hold dictionaries for each contig, storing length, gc and cov (to calculate weighted av. of gc and cov later for each cluster)
    cluster_contigs = dict() # Stores the cluster for each contig
    markers_in_cluster = dict() # Keyed by cluster and then PFAM, stores number of the PFAM in each cluster

    for i,row in master_table.iterrows():
        contig = row['contig']
        length = int(row['length'])
        cov = float(row['cov'])
        gc = float(row['gc'])
        cluster = row[args.column]

        contig_info[contig] = { 'length': length, 'cov': cov, 'gc': gc }
        cluster_contigs[contig] = cluster

        if cluster not in markers_in_cluster:
            markers_in_cluster[cluster] = dict()

        # Protect for instances where single_copy_PFAMs is empty, and gets converted to nan (a float) by pandas
        if isinstance(row['single_copy_PFAMs'], float):
            continue

        pfam_list = row['single_copy_PFAMs'].split(',')

        for pfam in pfam_list:
            if pfam not in markers_in_cluster[cluster]:
                markers_in_cluster[cluster][pfam] = 1
            else:
                markers_in_cluster[cluster][pfam] += 1

    # Load fasta file using biopython, and split into clusters
    cluster_sequences = dict() # Keyed by cluster, will hold lists of seq objects
    for seq_record in SeqIO.parse(args.fasta, 'fasta'):
        seq_name = str(seq_record.id)
        if seq_name in cluster_contigs:
            cluster = cluster_contigs[seq_name]
        else:
            continue

        if cluster not in cluster_sequences:
            cluster_sequences[cluster] = list()
        cluster_sequences[cluster].append(seq_record)

    # Output summary table plus individual fasta files
    summary_outfpath = os.path.join(args.outdir,'cluster_summary.tsv')
    fh = open(summary_outfpath, 'w')
    cols = [
        'cluster',
        'size',
        'longest_contig',
        'n50',
        'number_contigs',
        'completeness',
        'purity',
        'av_cov',
        'av_gc']
    header = '\t'.join(cols)+'\n'
    fh.write(header)

    for cluster in cluster_sequences:
        attributes = assess_assembly(cluster_sequences[cluster])
        total_size = attributes['size']
        longest_contig = attributes['largest_sequence']
        n50 = attributes['n50']
        number_contigs = attributes['number_sequences']

        if args.kingdom == 'bacteria':
            total_markers = 139
        elif args.kingdom == 'archaea':
            total_markers = 162

        number_markers_found = 0
        number_unique_markers = len(markers_in_cluster[cluster])
        for pfam in markers_in_cluster[cluster]:
            number_markers_found += markers_in_cluster[cluster][pfam]

        # The following protects for the edge case where a cluster has zero marker genes
        if number_unique_markers == 0:
            completness = 'unknown'
        else:
            completeness = (number_unique_markers / total_markers) * 100
        if number_markers_found == 0:
            purity = 'unknown'
        else:
            purity = (number_unique_markers / number_markers_found) * 100

        # Calculate average GC and cov, weighted by sequence length
        weighted_gc_av = 0.0
        weighted_cov_av = 0.0
        for seq_record in cluster_sequences[cluster]:
            seq_name = str(seq_record.id)
            seq_length = contig_info[seq_name]['length']
            seq_gc = contig_info[seq_name]['gc']
            seq_cov = contig_info[seq_name]['cov']
            seq_length_frac = seq_length / total_size
            weighted_gc_av += seq_gc * seq_length_frac
            weighted_cov_av += seq_cov * seq_length_frac

        # Write line in summary table
        line = '\t'.join(map(str,[
            cluster,
            total_size,
            longest_contig,
            n50,
            number_contigs,
            completeness,
            purity,
            weighted_cov_av,
            weighted_gc_av
        ]))+'\n'
        fh.write(line)

        # Write individual fasta file
        cluster_outfpath = os.path.join(args.outdir, f'cluster_{cluster}.fasta')
        SeqIO.write(cluster_sequences[cluster], cluster_outfpath, 'fasta')

    fh.close()

    # If the user has specified --do_taxonomy, then they also need to specify --db_dir
    if args.do_taxonomy:
        if not args.db_dir:
            logger.error('Error! If you want to analyze taxonomy, you need to specify a path to database files (--db_dir)')
            exit(1)

        if not os.path.isdir(args.db_dir):
            logger.error('Error! DB dir ' + args.db_dir + ' does not exist')
            exit(1)
        if not os.path.isfile(args.db_dir + '/names.dmp'):
            logger.error('Error! Cannot find names.dmp in ' + args.db_dir)
            exit(1)

        if not os.path.isfile(args.db_dir + '/nodes.dmp'):
            logger.error('Error! Cannot find nodes.dmp in ' + args.db_dir)
            exit(1)
        # Now run cluster_taxonomy.py
        taxonomy_output_path = os.path.join(args.outdir,'cluster_taxonomy.tab')
        cmd = f'cluster_taxonomy.py -t {args.bin_table} -c {args.column} -x {args.db_dir} -o {taxonomy_output_path}'

        run_command(cmd)


def main(args):
    logger.info(args.hello_world)
    # operations on args.positional
    # operations on args.optional

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='[%(asctime)s %(levelname)s] %(name)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    curdir = os.path.join(os.path.abspath(os.curdir), 'clusters')
    description = 'Script to summarize the assembly characteristics and taxonomy\
    of binned clusters, as well producing individual cluster fasta files'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-b', '--bin_table', metavar='<bin.tab>', help='path to the output from either run_autometa.py or ML_recruitment.py', required=True)
    parser.add_argument('-c', '--column', metavar='<bin column name>', help='the name of the column to use for binning purposes', default='cluster')
    parser.add_argument('-f', '--fasta', metavar='<assembly.fasta>', help='path to the assembly used to make the bin table', required=True)
    parser.add_argument('--outdir', help='</path/to/output/directory>', default=curdir)
    parser.add_argument('-k', '--kingdom', metavar='<archaea|bacteria>', help='kingdom to consider', choices=['bacteria', 'archaea'], default='bacteria')
    parser.add_argument('-t', '--do_taxonomy', help='carry out taxonomic analysis on the clusters (you must have already run make_taxonomy_table.py)', action='store_true')
    parser.add_argument('-db', '--db_dir', metavar='<dir>', help='Path to directory with taxdump files')
    args = vars(parser.parse_args())
    parser.add_argument('positional',help='<help text of positional arg>')
    parser.add_argument('--optional',help='<help text of optional arg>')
    parser.add_argument(
        '--hello-world',
        help='<help text of hello world>',
        default='Hello World')
    args = parser.parse_args()
    main(args)
