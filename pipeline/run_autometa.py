#!/usr/bin/env python

# Copyright 2018 Ian J. Miller, Evan Rees, Izaak Miller, Jason C. Kwan
#
# This file is part of Autometa.
#
# Autometa is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Autometa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with Autometa. If not, see <http://www.gnu.org/licenses/>.

import sys
import subprocess
import time
import logging
import os
import shutil

from argparse import ArgumentParser
from argparse import ArgumentDefaultsHelpFormatter


def init_logger(autom_path, db_path, out_path):
    # logger
    logger = logging.getLogger(out_path + "/run_autometa.py")
    hdlr = logging.FileHandler(out_path + "/run_autometa.log")
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    hdlr.setFormatter(formatter)
    logger.addHandler(hdlr)
    logger.setLevel(logging.DEBUG)

    # Output current branch and commit
    branch_command = "git -C " + autom_path + " branch | grep \* | sed 's/^..//'"
    branch = (
        subprocess.Popen(branch_command, shell=True, stdout=subprocess.PIPE)
        .communicate()[0]
        .rstrip()
    )

    commit_command = "git -C " + autom_path + " rev-parse --short HEAD"
    commit = (
        subprocess.Popen(commit_command, shell=True, stdout=subprocess.PIPE)
        .communicate()[0]
        .rstrip()
    )
    logger.info("Autometa branch: {}".format(branch))
    logger.info("Autometa commit: {}".format(commit))
    # Check programs
    progs = ["md5sum", "gzip", "tar", "gunzip", "wget"]
    for prog in progs:
        prog_version = subprocess.check_output([prog, "--version"]).split("\n")[0]
        logger.info("{}".format(prog_version))
    # Check 3rd party dependencies
    dmnd_v = subprocess.check_output(["diamond", "version"]).strip()
    logger.info("{}".format(dmnd_v))
    hmmscan_v = (
        subprocess.check_output(["hmmpress", "-h"]).split("\n")[1].replace("# ", "")
    )
    logger.info("{}".format(hmmscan_v))
    prodigal_v = (
        subprocess.Popen("prodigal -v", stderr=subprocess.PIPE, shell=True)
        .communicate()[1]
        .replace("\n", "")
    )
    logger.info("{}".format(prodigal_v))
    logger.info("DB Dir: {}".format(db_path))
    db_fpaths = os.listdir(db_path)
    for fpath in db_fpaths:
        logger.info(
            "DB (fname, size): {} {}".format(
                fpath, os.stat(db_path + "/" + fpath).st_size
            )
        )
    return logger


def run_command(command_string, stdout_path=None):
    # Function that checks if a command ran properly. If it didn't, then print an error message then quit
    logger.info("run_autometa.py, run_command: " + command_string)
    if stdout_path:
        f = open(stdout_path, "w")
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print("run_autometa.py: Error, the command:")
        print(command_string)
        print(("failed, with exit code " + str(exit_code)))
        exit(1)


def run_command_quiet(command_string):
    exit_code = subprocess.call(
        command_string, shell=True, stdout=FNULL, stderr=subprocess.STDOUT
    )

    if exit_code != 0:
        print("run_autometa.py: Error, the command:")
        print(command_string)
        print(("failed, with exit code " + str(exit_code)))
        exit(1)


def cythonize_lca_functions():
    logger.info(
        "{}/lca_functions.so not found, cythonizing lca_function.pyx for make_taxonomy_table.py".format(
            pipeline_path
        )
    )
    current_dir = os.getcwd()
    os.chdir(pipeline_path)
    run_command("python setup_lca_functions.py build_ext --inplace")
    os.chdir(current_dir)


def lca_compilation_check():
    lca_funcs_so = os.path.join(pipeline_path, "lca_functions.so")
    lca_fp = os.path.join(pipeline_path, "lca.py")
    if not os.path.isfile(lca_funcs_so):
        cythonize_lca_functions()
    elif os.path.getmtime(lca_funcs_so) < os.path.getmtime(lca_fp):
        print("lca.py updated. Recompiling lca functions.")
        build_dir = os.path.join(pipeline_path, "build")
        lca_funcs_c = lca_funcs_so.replace(".so", ".c")
        os.remove(lca_funcs_c)
        os.remove(lca_funcs_so)
        shutil.rmtree(build_dir)
        cythonize_lca_functions()
    else:
        print("lca_functions up-to-date")


def run_make_taxonomy_tab(fasta, length_cutoff):
    """Runs make_taxonomy_table.py and directs output to taxonomy.tab for run_autometa.py"""
    # Note we don't have to supply the cov_table here because earlier in this script we already run make_contig_table.py
    output_path = output_dir + "/taxonomy.tab"
    if cov_table:
        run_command(
            "{}/make_taxonomy_table.py -a {} -db {} -p {} -l {} -o {} -v {}".format(
                pipeline_path,
                fasta,
                db_dir_path,
                processors,
                length_cutoff,
                output_dir,
                cov_table,
            )
        )
    else:
        run_command(
            "{}/make_taxonomy_table.py -a {} -db {} -p {} -l {} -o {}".format(
                pipeline_path, fasta, db_dir_path, processors, length_cutoff, output_dir
            )
        )
    return output_path


def length_trim(fasta, length_cutoff):
    # will need to update path of this perl script
    outfile_name, ext = os.path.splitext(os.path.basename(fasta))
    outfile_name += ".filtered" + ext
    output_path = output_dir + "/" + outfile_name
    run_command(
        "{}/fasta_length_trim.pl {} {} {}".format(
            pipeline_path, fasta, length_cutoff, output_path
        )
    )
    return output_path


def make_contig_table(fasta, coverage_table=None):
    # Fasta is an absolute path
    output_filename, _ = os.path.splitext(os.path.basename(fasta))
    output_filename += ".tab"
    output_path = output_dir + "/" + output_filename
    if coverage_table:
        run_command(
            "{}/make_contig_table.py -a {} -c {} -o {}".format(
                pipeline_path, fasta, coverage_table, output_path
            )
        )
    else:
        run_command(
            "{}/make_contig_table.py -a {} -o {}".format(
                pipeline_path, fasta, output_path
            )
        )
    return output_path


def make_marker_table(fasta):
    if kingdom == "bacteria":
        hmm_marker_path = (
            autometa_path + "/single-copy_markers/Bacteria_single_copy.hmm"
        )
        hmm_cutoffs_path = (
            autometa_path + "/single-copy_markers/Bacteria_single_copy_cutoffs.txt"
        )
    elif kingdom == "archaea":
        hmm_marker_path = autometa_path + "/single-copy_markers/Archaea_single_copy.hmm"
        hmm_cutoffs_path = (
            autometa_path + "/single-copy_markers/Archaea_single_copy_cutoffs.txt"
        )

    # need to add processors to this script
    output_fname, _ = os.path.splitext(os.path.basename(fasta))
    output_fname += ".marker.tab"
    output_path = output_dir + "/" + output_fname
    if os.path.isfile(output_path):
        print("{} file already exists!".format(output_path))
        print("Continuing to next step...")
        logger.info("{} file already exists!".format(output_path))
        logger.info("Continuing to next step...")
    else:
        print(
            "Making the marker table with prodigal and hmmscan. This could take a while..."
        )
        logger.info(
            "Making the marker table with prodigal and hmmscan. This could take a while..."
        )
        if not os.path.exists((hmm_marker_path + ".h3m")):
            run_command_quiet("hmmpress -f {}".format(hmm_marker_path))
        run_command_quiet(
            "{}/make_marker_table.py -a {} -m {} -c {} -o {} -p {}".format(
                pipeline_path,
                fasta,
                hmm_marker_path,
                hmm_cutoffs_path,
                output_path,
                processors,
            )
        )
    return output_path


def recursive_dbscan(input_table, filtered_assembly, domain):
    recursive_dbscan_output_path = output_dir + "/recursive_dbscan_output.tab"
    k_mer_file = output_dir + "/k-mer_matrix"
    run_command(
        "{}/recursive_dbscan.py -t {} -a {} -d {} -k {}".format(
            pipeline_path, input_table, filtered_assembly, output_dir, domain
        )
    )

    return recursive_dbscan_output_path, k_mer_file


def combine_tables(table1_path, table2_path):
    comb_table_path = output_dir + "/combined_contig_info.tab"
    # Note: in this sub we assume that the tables have column 1 in common
    # Store lines of table 2, keyed by the value of the first column

    # First make data structure from table 2
    table2_lines = dict()
    with open(table2_path) as table2:
        for i, line in enumerate(table2):
            line_list = line.rstrip().split("\t")
            contig = line_list.pop(0)
            if i == 0:
                table_2_header = "\t".join(line_list)
            else:
                table2_lines[contig] = "\t".join(line_list)

    comb_table = open(comb_table_path, "w")
    with open(table1_path) as table1:
        for i, line in enumerate(table1):
            line_list = line.rstrip().split("\t")
            contig = line_list.pop(0)
            if i == 0:
                new_header = line.rstrip() + "\t" + table_2_header + "\n"
                comb_table.write(new_header)
            else:
                # We have to check whether the line exists in table 2
                if contig in table2_lines:
                    new_line = line.rstrip() + "\t" + table2_lines[contig] + "\n"
                    comb_table.write(new_line)

    comb_table.close()

    return comb_table_path


def ML_recruitment(input_table, matrix):
    ML_recruitment_output_path = output_dir + "/ML_recruitment_output.tab"
    run_command(
        "{}/ML_recruitment.py -t {} -p {} -r -m {} -o {}".format(
            pipeline_path, input_table, processors, matrix, ML_recruitment_output_path
        )
    )

    return ML_recruitment_output_path


pipeline_path = sys.path[0]
autometa_path = os.path.dirname(pipeline_path)

# argument parser
parser = ArgumentParser(
    description="Script to run the Autometa pipeline.",
    epilog="Please do not forget to cite us. Thank you for using Autometa!",
    formatter_class=ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    "-a",
    "--assembly",
    metavar="<assembly.fasta>",
    help="Path to metagenomic assembly fasta",
    required=True,
)
parser.add_argument(
    "-p",
    "--processors",
    metavar="<int>",
    help="Number of processors to use",
    type=int,
    default=1,
)
parser.add_argument(
    "-l",
    "--length_cutoff",
    metavar="<int>",
    help="Contig length cutoff to consider for binning in bp",
    default=10000,
    type=int,
)
parser.add_argument(
    "-c",
    "--completeness_cutoff",
    metavar="<float>",
    help="Completeness cutoff (in percent) to use for accepting clusters",
    type=float,
    default=20.0,
)
parser.add_argument(
    "-k",
    "--kingdom",
    metavar="<archaea|bacteria>",
    help="Kingdom to consider",
    choices=["bacteria", "archaea"],
    default="bacteria",
)
parser.add_argument(
    "-t",
    "--taxonomy_table",
    metavar="<taxonomy.tab>",
    help="Path to output of make_taxonomy_table.py",
)
parser.add_argument(
    "-o",
    "--output_dir",
    metavar="<dir>",
    help="Path to directory to store all output files",
    default=".",
)
parser.add_argument(
    "-r",
    "--ML_recruitment",
    help="Use ML to further recruit unclassified contigs",
    action="store_true",
)
parser.add_argument(
    "-m",
    "--maketaxtable",
    action="store_true",
    help="runs make_taxonomy_table.py before performing autometa binning. Must specify databases directory (-db)",
)
parser.add_argument(
    "-db",
    "--db_dir",
    metavar="<dir>",
    help="Path to directory with taxdump files. If this doesn't exist, the files will be automatically downloaded",
    required=False,
    default=autometa_path + "/databases",
)
parser.add_argument(
    "-v",
    "--cov_table",
    metavar="<coverage.tab>",
    help="Path to coverage table made by calculate_read_coverage.py. If this is not specified then coverage information will be extracted from contig names (SPAdes format)",
    required=False,
)

args = vars(parser.parse_args())

length_cutoff = args["length_cutoff"]
fasta_assembly = os.path.abspath(args["assembly"])
processors = args["processors"]
cluster_completeness = args["completeness_cutoff"]
kingdom = args["kingdom"].lower()
taxonomy_table_path = args["taxonomy_table"]
output_dir = os.path.abspath(args["output_dir"])
do_ML_recruitment = args["ML_recruitment"]
make_tax_table = args["maketaxtable"]
db_dir_path = os.path.abspath(args["db_dir"])
cov_table = args["cov_table"]

# Make output directory if it doesn't exist
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

logger = init_logger(autometa_path, db_dir_path, output_dir)

# check if fasta in path
if not os.path.isfile(fasta_assembly):
    print("Could not find {}...".format(fasta_assembly))
    logger.debug("Could not find {}...".format(fasta_assembly))
    exit(1)

# if make_tax_table specified but taxonomy_table_path not defined
if make_tax_table and not taxonomy_table_path:
    taxonomy_table_path = output_dir + "/taxonomy.tab"

# If coverage table is given, it must exist
if cov_table and not os.path.isfile(cov_table):
    print(("Error! Could not find coverage table at the following path: " + cov_table))
    exit(1)

# what input variables were and when you ran it (report fill path based on argparse)
logger.info("Input command: " + " ".join(sys.argv))

start_time = time.time()
FNULL = open(os.devnull, "w")

# run length trim and store output name
filtered_assembly = length_trim(fasta_assembly, length_cutoff)
if cov_table:
    contig_table = make_contig_table(filtered_assembly, cov_table)
else:
    contig_table = make_contig_table(filtered_assembly)
marker_tab_path = make_marker_table(filtered_assembly)

# Ensure lca functions are compiled and up-to-date
lca_compilation_check()

# Make combined table
if taxonomy_table_path and not make_tax_table:
    combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
elif taxonomy_table_path and make_tax_table:
    if not os.path.isfile(taxonomy_table_path):
        print(
            "Could not find {}, running make_taxonomy_table.py".format(
                taxonomy_table_path
            )
        )
        logger.debug(
            "Could not find {}, running make_taxonomy_table.py".format(
                taxonomy_table_path
            )
        )
        taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
        combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
    elif (
        os.path.isfile(taxonomy_table_path)
        and os.stat(taxonomy_table_path).st_size == 0
    ):
        print(
            "{} file is empty, running make_taxonomy_table.py".format(
                taxonomy_table_path
            )
        )
        logger.debug(
            "{} file is empty, running make_taxonomy_table.py".format(
                taxonomy_table_path
            )
        )
        taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
        combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
    else:
        print(
            "{} already exists, not performing make_taxonomy_table.py".format(
                taxonomy_table_path
            )
        )
        combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
elif not taxonomy_table_path and make_tax_table:
    taxonomy_table_path = run_make_taxonomy_tab(fasta_assembly, length_cutoff)
    combined_table_path = combine_tables(taxonomy_table_path, marker_tab_path)
else:
    combined_table_path = combine_tables(contig_table, marker_tab_path)

# If we carried out make_taxonomy_table.py, we now should use the appropriate kingdom bin instead of the raw
# input fasta
if make_tax_table:
    # First, check that the expected kingdom bin is there
    expected_kingdom_bin_path = output_dir + "/" + kingdom.title() + ".fasta"
    if (not os.path.isfile(expected_kingdom_bin_path)) or os.stat(
        expected_kingdom_bin_path
    ).st_size == 0:
        print(
            (
                "Error at make_taxonomy_table.py stage - file {0} is either not there or empty".format(
                    expected_kingdom_bin_path
                )
            )
        )
        exit(1)

    # Now change the input fasta to the output of make_taxonomy_table.py
    filtered_assembly = expected_kingdom_bin_path

recursive_dbscan_output, matrix_file = recursive_dbscan(
    combined_table_path, filtered_assembly, kingdom
)

if do_ML_recruitment:
    ML_recruitment(recursive_dbscan_output, matrix_file)

elapsed_time = time.strftime(
    "%H:%M:%S", time.gmtime(round((time.time() - start_time), 2))
)

print("Done!")
print("Elapsed time is {} (HH:MM:SS)".format(elapsed_time))
logger.info("Done!")
logger.info("Elapsed time is {} (HH:MM:SS)".format(elapsed_time))
FNULL.close()
