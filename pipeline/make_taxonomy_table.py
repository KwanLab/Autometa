#!/usr/bin/env python

# Copyright 2018 Ian J. Miller, Evan R. Rees, Izaak Miller, Jason C. Kwan
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

try:
    # python 2 compatible
    from urllib2 import urlopen
except ImportError:
    # python 3 compatible
    # https://stackoverflow.com/a/14510349/13118765
    from urllib.request import urlopen as urllib_urlopen
    from urllib.request import Request

    def urlopen(url):
        return urllib_urlopen(Request(url))


import subprocess
import os
import shutil

import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO


PIPELINE = os.path.dirname(os.path.realpath(__file__))
AUTOMETA_DATABASES = os.path.join(os.path.dirname(PIPELINE), "databases")


def run_command(command_string, stdout_path=None):
    # Function that checks if a command ran properly
    # If it didn't, then print an error message then quit
    print("make_taxonomy_table.py, run_command: " + command_string)
    if stdout_path:
        f = open(stdout_path, "w")
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print("make_taxonomy_table.py: Error, the command:")
        print(command_string)
        print("failed, with exit code " + str(exit_code))
        exit(1)


def run_command_return(command_string, stdout_path=None):
    # Function that checks if a command ran properly.
    # If it didn't, then print an error message then quit
    print("make_taxonomy_table.py, run_command: " + command_string)
    if stdout_path:
        f = open(stdout_path, "w")
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    return exit_code


def cythonize_lca_functions():
    lca_functions_so = os.path.join(PIPELINE, "lca_functions.so")
    print(
        "{} not found, cythonizing lca_function.pyx for make_taxonomy_table.py".format(
            lca_functions_so
        )
    )
    current_dir = os.getcwd()
    os.chdir(PIPELINE)
    run_command("python setup_lca_functions.py build_ext --inplace")
    os.chdir(current_dir)


def lca_compilation_check():
    lca_funcs_so = os.path.join(PIPELINE, "lca_functions.so")
    lca_fp = os.path.join(PIPELINE, "lca.py")
    if not os.path.isfile(lca_funcs_so):
        cythonize_lca_functions()
    elif os.path.getmtime(lca_funcs_so) < os.path.getmtime(lca_fp):
        print("lca.py updated. Recompiling lca functions.")
        build_dir = os.path.join(PIPELINE, "build")
        lca_funcs_c = lca_funcs_so.replace(".so", ".c")
        os.remove(lca_funcs_c)
        os.remove(lca_funcs_so)
        shutil.rmtree(build_dir)
        cythonize_lca_functions()
    else:
        print("lca_functions up-to-date")


def download_file(destination_dir, file_url, md5_url):
    filename = os.path.basename(file_url)
    md5name = os.path.basename(md5_url)

    while True:
        outfpath = os.path.join(destination_dir, filename)
        md5_fpath = os.path.join(destination_dir, md5name)
        run_command("wget {} -O {}".format(file_url, outfpath))
        run_command("wget {} -O {}".format(md5_url, md5_fpath))

        downloaded_md5 = subprocess.check_output(["md5sum", outfpath]).split(" ")[0]

        with open(md5_fpath, "r") as check_md5_file:
            check_md5 = check_md5_file.readline().split(" ")[0]

        if downloaded_md5 == check_md5:
            print("md5 checksum successful. Continuing...")
            break
        else:
            print("md5 checksum unsuccessful. Retrying...")


def md5IsCurrent(local_md5_path, remote_md5_url):
    remote_md5_handle = urlopen(remote_md5_url)
    # python3 request returns binary encoded file. Therefore we need to explicitly decode to utf-8
    remote_md5 = remote_md5_handle.readline().decode("utf-8").split(" ")[0]
    with open(local_md5_path, "r") as local_md5_file:
        local_md5 = local_md5_file.readline().split(" ")[0]
    return local_md5 == remote_md5


def prepare_databases(outdir, db="all", update=False):
    """Updates databases for AutoMeta usage"""
    # Downloading files for db population
    if db == "all" or db == "nr":
        nr_db_url = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz"
        nr_db_md5_url = nr_db_url + ".md5"
        # First download nr if we don't yet have it OR it is not up to date
        nr_md5_fpath = os.path.join(outdir, "nr.gz.md5")
        if os.path.isfile(nr_md5_fpath) and os.path.getsize(nr_md5_fpath):
            # Since a checksum exists, we should determine whether the nr database is up-to-date.
            # Only update if the update argument is True.
            # Otherwise, notify the user when the database is not up-to-date.
            if update and not md5IsCurrent(nr_md5_fpath, nr_db_md5_url):
                print("md5 is not current. Updating nr.dmnd")
                download_file(outdir, nr_db_url, nr_db_md5_url)
            elif not update and not md5IsCurrent(nr_md5_fpath, nr_db_md5_url):
                print(
                    "md5 is not current. Consider supplying the --update argument to update the nr.dmnd database"
                )
            else:
                print("nr.dmnd database is up to date.")
        else:
            print("{} not found. Downloading nr.gz".format(nr_md5_fpath))
            download_file(outdir, nr_db_url, nr_db_md5_url)

        nr_fpath = os.path.join(outdir, "nr.gz")
        # Now we make the diamond database
        print("building nr.dmnd database, this may take some time")
        cmd = " ".join(
            [
                "diamond makedb",
                "--in {}".format(nr_fpath),
                "--db {}".format(nr_fpath.rstrip(".gz")),
                "-p {}".format(num_processors),
            ]
        )
        returnCode = subprocess.call(cmd, shell=True)
        if returnCode == 0:  # i.e. job was successful
            # Make an md5 file to signal that we have built the database successfully
            os.remove(nr_fpath)
            print("nr.dmnd updated")
        else:
            print("nr.dmnd update FAILED!")
    if db == "all" or db == "acc2taxid":
        accession2taxid_url = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
        accession2taxid_md5_url = accession2taxid_url + ".md5"
        # Download prot.accession2taxid.gz only if the version we have is not current
        acc2taxid_md5_fpath = os.path.join(outdir, "prot.accession2taxid.gz.md5")
        if os.path.isfile(acc2taxid_md5_fpath) and os.path.getsize(acc2taxid_md5_fpath):
            if update and not md5IsCurrent(
                acc2taxid_md5_fpath, accession2taxid_md5_url
            ):
                print("md5 is not current. Updating prot.accession2taxid")
                download_file(outdir, accession2taxid_url, accession2taxid_md5_url)
            if not update and not md5IsCurrent(
                acc2taxid_md5_fpath, accession2taxid_md5_url
            ):
                print(
                    "md5 is not current. Consider updating prot.accession2taxid with the --update flag."
                )
            else:
                print("md5 is up-to-date for prot.accession2taxid")
        else:
            print("updating prot.accession2taxid")
            download_file(outdir, accession2taxid_url, accession2taxid_md5_url)

        acc2taxid_fpath = os.path.join(outdir, "prot.accession2taxid.gz")
        if os.path.isfile(acc2taxid_fpath):
            print(
                "Gunzipping prot.accession2taxid gzipped file\nThis may take some time..."
            )
            cmd = "gunzip -9vNf {}".format(acc2taxid_fpath)
            run_command(cmd, stdout_path=acc2taxid_fpath.rstrip(".gz"))
            print("prot.accession2taxid updated")
    if db == "all" or db == "taxdump":
        taxdump_url = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
        taxdump_md5_url = taxdump_url + ".md5"
        # Download taxdump only if the version we have is not current
        taxdump_md5_fpath = os.path.join(outdir, "taxdump.tar.gz.md5")

        if os.path.isfile(taxdump_md5_fpath) and os.path.getsize(taxdump_md5_fpath):
            if update and not md5IsCurrent(taxdump_md5_fpath, taxdump_md5_url):
                print("updating nodes.dmp, names.dmp, merged.dmp")
                download_file(outdir, taxdump_url, taxdump_md5_url)
            if not update and not md5IsCurrent(taxdump_md5_fpath, taxdump_md5_url):
                print(
                    "taxdump files are not current, considering updating with --update flag."
                )
            else:
                print("taxdump files are up-to-date.")
        else:
            print("updating nodes.dmp, names.dmp, merged.dmp")
            download_file(outdir, taxdump_url, taxdump_md5_url)

        taxdump_fpath = os.path.join(outdir, "taxdump.tar.gz")
        if os.path.isfile(taxdump_fpath):
            run_command(
                "tar -xzf {} -C {} names.dmp nodes.dmp merged.dmp".format(
                    taxdump_fpath, outdir
                )
            )
            os.remove(taxdump_fpath)
            print("nodes.dmp, names.dmp, merged.dmp updated")


def check_dbs(db_path, update=False):
    """
    Determines what files need to be downloaded/updated depending on
    what database path was specified
    """
    if os.path.realpath(db_path) == AUTOMETA_DATABASES:
        db_dict = {
            "nr": ["nr.dmnd", "nr.gz.md5"],
            "acc2taxid": ["prot.accession2taxid.gz.md5", "prot.accession2taxid"],
            "taxdump": ["merged.dmp", "names.dmp", "nodes.dmp", "taxdump.tar.gz.md5"],
        }
    else:
        db_dict = {
            "nr": ["nr.dmnd"],
            "acc2taxid": ["prot.accession2taxid"],
            "taxdump": ["names.dmp", "nodes.dmp", "merged.dmp"],
        }
    db_files = os.listdir(db_path)
    for db, fpaths in db_dict.items():
        for fpath in fpaths:
            if fpath not in db_files:
                print(
                    "{0} database not found, downloading/formatting.\n\
				This may take some time...".format(
                        db
                    )
                )
                prepare_databases(outdir=db_path, db=db, update=update)


def length_trim(fasta_path, length_cutoff):
    input_fname, ext = os.path.splitext(os.path.basename(fasta_path))
    # Trim the length of fasta file
    outfname = input_fname + ".filtered" + ext
    outfile_path = os.path.join(output_dir, outfname)
    script = os.path.join(PIPELINE, "fasta_length_trim.pl")
    cmd = " ".join(map(str, [script, fasta_path, length_cutoff, outfile_path]))
    run_command(cmd)
    return outfile_path


def run_prodigal(path_to_assembly):
    assembly_fname, _ = os.path.splitext(os.path.basename(path_to_assembly))
    output_path = os.path.join(output_dir, assembly_fname + ".orfs.faa")
    if os.path.isfile(output_path):
        print("{} file already exists!".format(output_path))
        print("Continuing to next step...")
    else:
        run_command(
            "prodigal -i {} -a {}/{}.orfs.faa -p meta -m -o {}/{}.txt".format(
                path_to_assembly, output_dir, assembly_fname, output_dir, assembly_fname
            )
        )


def run_diamond(orfs_fpath, diamond_db_path, num_processors, outfpath):
    view_output = orfs_fpath + ".blastp"
    tmp_dir_path = os.path.join(os.path.dirname(orfs_fpath), "tmp")
    if not os.path.isdir(tmp_dir_path):
        os.makedirs(
            tmp_dir_path
        )  # This will give an error if the path exists but is a file instead of a dir
    cmds = [
        "diamond blastp",
        "--evalue 1e-5",
        "--max-target-seqs 200",
        "--outfmt 6",
        "--query {}.faa".format(orfs_fpath),
        "--db {}".format(diamond_db_path),
        "-p {}".format(num_processors),
        "--out {}".format(outfpath),
        "-t {}".format(tmp_dir_path),
    ]
    cmd = " ".join(cmds)
    error = run_command_return(cmd)
    # If there is an error, attempt to rebuild NR
    if error == 134 or error == str(134):
        print("Fatal: Not enough disk space for diamond alignment archive!")
        exit(1)
    if error:
        print(
            "Error when performing diamond blastp:\n{}\nAttempting to correct by rebuilding nr...".format(
                error
            )
        )
        prepare_databases(outdir=db_dir_path, db="nr", update=False)
        # Retry with rebuilt nr.dmnd
        run_command(cmd)

    return outfpath


# blast2lca using accession numbers#
def run_blast2lca(input_file, taxdump_path):
    fname = os.path.splitext(os.path.basename(input_file))[0] + ".lca"
    output = os.path.join(output_dir, fname)
    if os.path.isfile(output) and os.path.getsize(output):
        print("{} file already exists! Continuing to next step...".format(output))
    else:
        lca_script = os.path.join(PIPELINE, "lca.py")
        cmd = "{} database_directory {} {}".format(lca_script, db_dir_path, input_file)
        run_command(cmd)
    return output


def run_taxonomy(
    assembly_path,
    lca_fpath,
    db_dir_path,
    coverage_table,
    bgcs_path=None,
    orfs_path=None,
):  # Have to update this
    assembly_fname, _ = os.path.splitext(os.path.basename(assembly_path))
    contig_tab_fpath = os.path.join(output_dir, assembly_fname + ".tab")
    contig_table_script = os.path.join(PIPELINE, "make_contig_table.py")
    # Only make the contig table if it doesn't already exist
    if not os.path.isfile(contig_tab_fpath):
        cmd = "{} -a {} -o {}".format(
            contig_table_script, assembly_path, contig_tab_fpath
        )
        if coverage_table:
            cmd += " -c {}".format(coverage_table)
        elif single_genome_mode:
            cmd += " -n"
        run_command(cmd)
    if bgcs_path:
        mask_bgcs_script = os.path.join(PIPELINE, "mask_bgcs.py2.7")
        cmd = "{} --bgc {} --orfs {} --lca {}"
        cmd = cmd.format(mask_bgcs_script, bgcs_path, orfs_path, lca_fpath)
        run_command(cmd)
        unmasked_fname = os.path.basename(lca_fpath).replace(".lca", ".unmasked.tsv")
        lca_fpath = os.path.join(output_dir, unmasked_fname)
    # two_files_generated: *.masked.tsv, *.unmasked.tsv
    add_contig_taxa_script = os.path.join(PIPELINE, "add_contig_taxonomy.py")
    taxonomy_fp = os.path.join(output_dir, "taxonomy.tab")
    cmd = " ".join(
        [add_contig_taxa_script, contig_tab_fpath, lca_fpath, db_dir_path, taxonomy_fp]
    )
    run_command(cmd)
    return taxonomy_fp


# argument parser
parser = ArgumentParser(
    description="Script to generate the contig taxonomy table.",
    epilog="Output will be directed to recursive_dbscan.py",
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
    help="Number of processors to use.",
    type=int,
    default=1,
)
parser.add_argument(
    "-db",
    "--db_dir",
    metavar="<dir>",
    help="Path to directory with taxdump, protein accessions and diamond (NR) \
	protein files. If this path does not exist, will create and download files.",
    required=False,
    default=AUTOMETA_DATABASES,
)
parser.add_argument(
    "-udb",
    "--user_prot_db",
    metavar="<user_prot_db>",
    help="Replaces the default diamond database (nr.dmnd)",
    required=False,
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
    "-v",
    "--cov_table",
    metavar="<coverage.tab>",
    help="Path to coverage table made by calculate_read_coverage.py. If this is \
	not specified then coverage information will be extracted from contig names (SPAdes format)",
    required=False,
)
parser.add_argument(
    "-o",
    "--output_dir",
    metavar="<dir>",
    help="Path to directory to store output files",
    default=".",
)
parser.add_argument(
    "-bgc",
    "--bgcs_dir",
    metavar="<dir>",
    help="Path to directory of biosynthetic gene clusters. Masks BGCs",
)
parser.add_argument(
    "-s", "--single_genome", help="Specifies single genome mode", action="store_true"
)
parser.add_argument(
    "-u",
    "--update",
    required=False,
    action="store_true",
    help="Checks/Adds/Updates: nodes.dmp, names.dmp, merged.dmp, accession2taxid, nr.dmnd files within specified directory.",
)

args = vars(parser.parse_args())

db_dir_path = os.path.abspath(args["db_dir"])
usr_prot_path = args["user_prot_db"]
num_processors = args["processors"]
length_cutoff = args["length_cutoff"]
fasta_path = args["assembly"]
cov_table = args["cov_table"]
output_dir = args["output_dir"]
single_genome_mode = args["single_genome"]

bgcs_dir = args["bgcs_dir"]
fasta_fname, _ = os.path.splitext(os.path.basename(fasta_path))
prodigal_output = os.path.join(output_dir, "{}.filtered.orfs".format(fasta_fname))
diamond_outfpath = prodigal_output + ".blastp"

# If cov_table defined, we need to check the file exists
if cov_table:
    if not os.path.isfile(cov_table):
        print(
            "Error! Could not find coverage table at the following path: " + cov_table
        )
        exit(1)

# Check that output dir exists, and create it if it doesn't
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

lca_compilation_check()

if not os.path.isdir(db_dir_path):
    # Verify the 'Autometa databases' directory exists
    print(
        "No databases directory found, creating and populating AutoMeta databases directory\n\
	This may take some time..."
    )
    os.mkdir(db_dir_path)
    prepare_databases(outdir=db_dir_path, db="all", update=args["update"])
elif not os.listdir(db_dir_path):
    # The 'Autometa databases' directory is empty
    print(
        "AutoMeta databases directory empty, populating with appropriate databases.\n\
	This may take some time..."
    )
    prepare_databases(outdir=db_dir_path, db="all", update=args["update"])
else:
    check_dbs(db_path=db_dir_path, update=args["update"])

names_dmp_path = os.path.join(db_dir_path, "names.dmp")
nodes_dmp_path = os.path.join(db_dir_path, "nodes.dmp")
accession2taxid_path = os.path.join(db_dir_path, "prot.accession2taxid")
diamond_db_path = os.path.join(db_dir_path, "nr.dmnd")
current_taxdump_md5 = os.path.join(db_dir_path, "taxdump.tar.gz.md5")
current_acc2taxid_md5 = os.path.join(db_dir_path, "prot.accession2taxid.gz.md5")
current_nr_md5 = os.path.join(db_dir_path, "nr.gz.md5")

if usr_prot_path:
    usr_prot_path = os.path.abspath(usr_prot_path)
    if os.path.isdir(usr_prot_path):
        print(
            "You have provided a directory {}. \
		--user_prot_db requires a file path.".format(
                usr_prot_path
            )
        )
        exit(1)
    elif not os.path.isfile(usr_prot_path):
        print("{} is not a file.".format(usr_prot_path))
        exit(1)
    else:
        diamond_db_path = usr_prot_path

if args["update"]:
    print("Checking database directory for updates")
    prepare_databases(outdir=db_dir_path, db="all", update=args["update"])

filtered_assembly = os.path.join(output_dir, "{}.filtered.fasta".format(fasta_fname))
if not os.path.isfile(filtered_assembly):
    filtered_assembly = length_trim(fasta_path, length_cutoff)

if not os.path.isfile(prodigal_output + ".faa"):
    print("Prodigal output not found. Running prodigal...")
    # Check for file and if it doesn't exist run make_marker_table
    run_prodigal(filtered_assembly)

if not os.path.isfile(diamond_outfpath):
    print("Could not find {}. Running diamond blast... ".format(diamond_outfpath))
    diamond_output = run_diamond(
        prodigal_output, diamond_db_path, num_processors, diamond_outfpath
    )
elif not os.path.getsize(diamond_outfpath):
    print("{} file is empty. Re-running diamond blast...".format(diamond_outfpath))
    diamond_output = run_diamond(
        prodigal_output, diamond_db_path, num_processors, diamond_outfpath
    )
elif not os.path.isfile(diamond_outfpath):
    print("{} not found. \nExiting...".format(diamond_outfpath))
    exit(1)
else:
    diamond_output = diamond_outfpath

if not os.path.isfile(prodigal_output + ".lca"):
    print("Could not find {}. Running lca...".format(prodigal_output + ".lca"))
    blast2lca_output = run_blast2lca(diamond_output, db_dir_path)
elif not os.path.getsize(prodigal_output + ".lca"):
    print("{} file is empty. Re-running lca...".format(prodigal_output + ".lca"))
    blast2lca_output = run_blast2lca(diamond_output, db_dir_path)
else:
    blast2lca_output = prodigal_output + ".lca"

taxonomy_table = os.path.join(output_dir, "taxonomy.tab")
if not os.path.isfile(taxonomy_table) or not os.path.getsize(taxonomy_table):
    print("Running add_contig_taxonomy.py... ")
    if bgcs_dir:
        taxonomy_table = run_taxonomy(
            assembly_path=filtered_assembly,
            lca_fpath=blast2lca_output,
            db_dir_path=db_dir_path,
            coverage_table=cov_table,
            bgcs_path=bgcs_dir,
            orfs_path=prodigal_output + ".faa",
        )
    else:
        taxonomy_table = run_taxonomy(
            assembly_path=filtered_assembly,
            lca_fpath=blast2lca_output,
            db_dir_path=db_dir_path,
            coverage_table=cov_table,
        )
else:
    print("taxonomy.tab exists... Splitting original contigs into kingdoms")

# Split the original contigs into sets for each kingdom
taxonomy_pd = pd.read_table(taxonomy_table)
categorized_seq_objects = {}
all_seq_records = {}

# Load fasta file
for seq_record in SeqIO.parse(filtered_assembly, "fasta"):
    all_seq_records[seq_record.id] = seq_record

for i, row in taxonomy_pd.iterrows():
    kingdom = row["kingdom"]
    contig = row["contig"]
    if contig not in all_seq_records:
        # Using filtered assembly, taxonomy.tab contains contigs not filtered
        print("{0} below length filter, skipping.".format(contig))
        continue
    if kingdom in categorized_seq_objects:
        categorized_seq_objects[kingdom].append(all_seq_records[contig])
    else:
        categorized_seq_objects[kingdom] = [all_seq_records[contig]]

# Now we write the component fasta files
if not single_genome_mode:
    for kingdom in categorized_seq_objects:
        seq_list = categorized_seq_objects[kingdom]
        output_path = os.path.join(output_dir, "{}.fasta".format(kingdom))
        SeqIO.write(seq_list, output_path, "fasta")

print("Done!")
