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
    from urllib.request import urlopen
except ImportError:
    # python 3 compatible
    # https://stackoverflow.com/a/14510349/13118765
    from urllib.request import urlopen as urllib_urlopen
    from urllib.request import Request

    def urlopen(url):
        return urllib_urlopen(Request(url))

import gzip
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
    print(("make_taxonomy_table.py, run_command: " + command_string))
    if stdout_path:
        f = open(stdout_path, "w")
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print("make_taxonomy_table.py: Error, the command:")
        print(command_string)
        print(("failed, with exit code " + str(exit_code)))
        exit(1)


def run_command_return(command_string, stdout_path=None):
    # Function that checks if a command ran properly.
    # If it didn't, then print an error message then quit
    print(("make_taxonomy_table.py, run_command: " + command_string))
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
        (
            "{} not found, cythonizing lca_function.pyx for make_taxonomy_table.py".format(
                lca_functions_so
            )
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
        run_command(f"wget {file_url} -O {outfpath}")
        run_command(f"wget {md5_url} -O {md5_fpath}")

        downloaded_md5 = subprocess.check_output(["md5sum", outfpath]).decode("utf-8").split(" ")[0]

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
        nr_db_md5_url = f"{nr_db_url}.md5"
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
            print(f"{nr_md5_fpath} not found. Downloading nr.gz")
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
            run_command(f"tar -xzf {taxdump_fpath} -C {outdir} names.dmp nodes.dmp merged.dmp")
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
            "acc2taxid": ["prot.accession2taxid.gz.md5", "prot.accession2taxid.gz"],
            "taxdump": ["merged.dmp", "names.dmp", "nodes.dmp", "taxdump.tar.gz.md5"],
        }
    else:
        db_dict = {
            "nr": ["nr.dmnd"],
            "acc2taxid": ["prot.accession2taxid.gz"],
            "taxdump": ["names.dmp", "nodes.dmp", "merged.dmp"],
        }
    db_files = os.listdir(db_path)
    for db, fpaths in list(db_dict.items()):
        for fpath in fpaths:
            if db == "acc2taxid":
                fpath = (
                    fpath.strip(".gz") if os.path.exists(fpath.strip(".gz")) else fpath
                )
            if fpath not in db_files:
                print(
                    f"{db} database not found, downloading/formatting... This may take some time..."
                )
                prepare_databases(outdir=db_path, db=db, update=update)


def length_trim(fasta, length_cutoff, output):
    filename = os.path.basename(fasta)
    filename = filename.strip(".gz") if filename.endswith(".gz") else filename
    infh = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta)
    records = [
        record
        for record in SeqIO.parse(infh, "fasta")
        if len(record.seq) >= length_cutoff
    ]
    infh.close()
    n_seqs = SeqIO.write(records, output, "fasta")
    print(f"Wrote length-filtered assembly ({length_cutoff}bp cutoff) to {output} ({n_seqs:,} seqs written)")
    return output


def run_prodigal(path_to_assembly, output):
    assembly_fname, _ = os.path.splitext(os.path.basename(path_to_assembly))
    if os.path.isfile(output):
        print(f"{output} file already exists!")
        print("Continuing to next step...")
    else:
        run_command(f"prodigal -i {path_to_assembly} -a {output} -p meta -m -o {output_dir}/{assembly_fname}.prodigal.txt")
        print(f"Wrote prodigal-called ORFs to {output}")
    return output

def run_diamond(orfs_fpath, diamond_db_path, num_processors, outfpath):
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
        f"--query {orfs_fpath}",
        f"--db {diamond_db_path}",
        f"-p {num_processors}",
        f"--out {outfpath}",
        f"-t {tmp_dir_path}",
    ]
    cmd = " ".join(cmds)
    error = run_command_return(cmd)
    # If there is an error, attempt to rebuild NR
    if error == 134 or error == str(134):
        print("Fatal: Not enough disk space for diamond alignment archive!")
        exit(1)
    if error:
        print(
            f"Error when performing diamond blastp:\n{error}\nAttempting to correct by rebuilding nr..."
        )
        prepare_databases(outdir=db_dir_path, db="nr", update=False)
        # Retry with rebuilt nr.dmnd
        run_command(cmd)

    return outfpath


# blast2lca using accession numbers#
def run_blast2lca(input_file, db_dir_path):
    fname = os.path.splitext(os.path.basename(input_file))[0] + ".lca"
    output = os.path.join(output_dir, fname)
    if os.path.isfile(output) and os.path.getsize(output):
        print(f"{output} file already exists! Continuing to next step...")
    else:
        lca_script = os.path.join(PIPELINE, "lca.py")
        cmd = f"{lca_script} database_directory {db_dir_path} {input_file}"
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
        cmd = f"{contig_table_script} -a {assembly_path} -o {contig_tab_fpath}"
        if coverage_table:
            cmd += f" -c {coverage_table}"
        elif single_genome_mode:
            cmd += " -n"
        run_command(cmd)
    if bgcs_path:
        mask_bgcs_script = os.path.join(PIPELINE, "mask_bgcs.py2.7")
        cmd = f"{mask_bgcs_script} --bgc {bgcs_path} --orfs {orfs_path} --lca {lca_fpath}"
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


# If cov_table defined, we need to check the file exists
if cov_table:
    if not os.path.isfile(cov_table):
        print(f"Error! Could not find coverage table at the following path: {cov_table}")
        exit(1)

# Check that output dir exists, and create it if it doesn't
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

lca_compilation_check()

if not os.path.isdir(db_dir_path):
    # Verify the 'Autometa databases' directory exists
    print("No databases directory found, creating and populating AutoMeta databases directory\nThis may take some time...")
    os.mkdir(db_dir_path)
    prepare_databases(outdir=db_dir_path, db="all", update=args["update"])
elif not os.listdir(db_dir_path):
    # The 'Autometa databases' directory is empty
    print("Autometa databases directory empty, populating with appropriate databases.\nThis may take some time...")
    prepare_databases(outdir=db_dir_path, db="all", update=args["update"])
else:
    check_dbs(db_path=db_dir_path, update=args["update"])

names_dmp_path = os.path.join(db_dir_path, "names.dmp")
nodes_dmp_path = os.path.join(db_dir_path, "nodes.dmp")
accession2taxid_path = os.path.join(db_dir_path, "prot.accession2taxid")
if not os.path.exists(accession2taxid_path) and os.path.exists(
    f"{accession2taxid_path}.gz"
):
    accession2taxid_path = f"{accession2taxid_path}.gz"
diamond_db_path = os.path.join(db_dir_path, "nr.dmnd")
current_taxdump_md5 = os.path.join(db_dir_path, "taxdump.tar.gz.md5")
current_acc2taxid_md5 = os.path.join(db_dir_path, "prot.accession2taxid.gz.md5")
current_nr_md5 = os.path.join(db_dir_path, "nr.gz.md5")

if usr_prot_path:
    usr_prot_path = os.path.abspath(usr_prot_path)
    if os.path.isdir(usr_prot_path):
        print(f"You have provided a directory {usr_prot_path}. --user_prot_db requires a file path.")
        exit(1)
    elif not os.path.isfile(usr_prot_path):
        print(f"{usr_prot_path} is not a file.")
        exit(1)
    else:
        diamond_db_path = usr_prot_path

if args["update"]:
    print("Checking database directory for updates")
    prepare_databases(outdir=db_dir_path, db="all", update=args["update"])

fasta_fname, _ = os.path.splitext(os.path.basename(fasta_path).replace(".gz", ""))
filtered_assembly = os.path.join(output_dir, f"{fasta_fname}.filtered.fna")
if not os.path.isfile(filtered_assembly):
    filtered_assembly = length_trim(fasta=fasta_path, length_cutoff=length_cutoff, output=filtered_assembly)

assembly_fname, ext = os.path.splitext(os.path.basename(filtered_assembly))
prodigal_output = os.path.join(output_dir, f"{assembly_fname}.orfs.faa")
diamond_outfpath = os.path.join(output_dir, f"{assembly_fname}.orfs.blastp")

if not os.path.isfile(prodigal_output):
    print("Prodigal output not found. Running prodigal...")
    # Check for file and if it doesn't exist run make_marker_table
    prodigal_output = run_prodigal(path_to_assembly=filtered_assembly, output=prodigal_output)

if not os.path.isfile(diamond_outfpath):
    print(f"Could not find {diamond_outfpath}. Running diamond blast... ")
    diamond_output = run_diamond(
        prodigal_output, diamond_db_path, num_processors, diamond_outfpath
    )
elif not os.path.getsize(diamond_outfpath):
    print(f"{diamond_outfpath} file is empty. Re-running diamond blast...")
    diamond_output = run_diamond(
        prodigal_output, diamond_db_path, num_processors, diamond_outfpath
    )
elif not os.path.isfile(diamond_outfpath):
    print(f"{diamond_outfpath} not found. \nExiting...")
    exit(1)
else:
    diamond_output = diamond_outfpath

lca_outfpath = os.path.join(output_dir, f"{assembly_fname}.orfs.lca")
if not os.path.isfile(lca_outfpath):
    print(f"Could not find {lca_outfpath}. Running lca...")
    blast2lca_output = run_blast2lca(diamond_output, db_dir_path)
elif not os.path.getsize(lca_outfpath):
    print(f"{lca_outfpath} file is empty. Re-running lca...")
    blast2lca_output = run_blast2lca(diamond_output, db_dir_path)
else:
    blast2lca_output = lca_outfpath

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
            orfs_path=prodigal_output,
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
        print(("{0} below length filter, skipping.".format(contig)))
        continue
    if kingdom in categorized_seq_objects:
        categorized_seq_objects[kingdom].append(all_seq_records[contig])
    else:
        categorized_seq_objects[kingdom] = [all_seq_records[contig]]

# Now we write the component fasta files
if not single_genome_mode:
    for kingdom in categorized_seq_objects:
        seq_list = categorized_seq_objects[kingdom]
        output_path = os.path.join(output_dir, f"{kingdom}.fasta")
        SeqIO.write(seq_list, output_path, "fasta")

print("Done!")
