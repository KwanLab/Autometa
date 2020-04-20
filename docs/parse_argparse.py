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

Script used to parse the arparse block of all the autometa scripts, copy them to a new
file and run the help text (--help) from there.
Used as an alternative to sphinxcontrib-programoutput which used all the heavy dependencies
preventing the intergration with readthedocs.
"""

import glob
import subprocess
import os
import textwrap
import shutil
import tempfile


def parse_files():
    """
    Makes a list of all the scripts that needs to be documented

    Returns
    -------
    list
        file paths that needs to be documented
    """
    list_fpath = []
    for fpath in glob.glob(os.path.join(os.path.dirname(os.getcwd()), "autometa", "**", "*.py"), recursive=True):
        if "__" in fpath:  # to not include __init__ files
            continue
        if os.path.basename(os.path.dirname(fpath)) == "validation":
            continue
        list_fpath.append(fpath)
    return(list_fpath)


def get_argparse_block(list_fpath):
    """
    Copies the argpasrse block of each script in a temporary file.

    Parameters
    ----------
    list_fpath : list
        path of scripts to be documented

    Returns
    -------
    list
        temporary files with the argparse block
    """
    list_path_temp = []
    for fpath in list_fpath:
        with open(fpath) as fh_parse:
            fd_temp, path_temp = tempfile.mkstemp()
            # path_temp = "/home/the_bio_informatician/Autometa/docs/temp_write.py"
            with open(path_temp, "w") as fh_write:
                writing = False
                for line in fh_parse:
                    line = line.strip()
                    if line == "import argparse":
                        fh_write.write(line + "\n" + "import os \n" +
                                       "import multiprocessing as mp \n" + "import work_queue as wq \n")
                        writing = True
                    # Add usage  = ScriptName.py in argparse block
                    if "parser = " in line:
                        fname = os.path.basename(fpath)
                        # Uses the starting `(` of Argument.parser, and replaces it with ( usage=
                        usage = f"( usage = \"{fname}\", "
                        line = line.replace("(", usage)
                    if writing:  # makes  DEFAULT_FPATH a `str` in user.py
                        if "DEFAULT_FPATH" in line:
                            line = line.replace(
                                "DEFAULT_FPATH", "\"DEFAULT_FPATH\"")
                    if writing:
                        # above if statment writes config."DEFAULT_FPATH" in database.py
                        # this if statement corrects it to "DEFAULT_FPATH"
                        if "config.\"DEFAULT_FPATH" in line:
                            line = line.replace(
                                "config.\"DEFAULT_FPATH", "\"config.DEFAULT_FPATH")
                    if writing:
                        fh_write.write(line + "\n")
                    if line == "args = parser.parse_args()":
                        writing = False
        list_path_temp.append(path_temp)
    return list_path_temp


def path_dir_docs_scripts(list_fpath):
    """
    Replaces the autometa directory in the path of scripts that needs to be documented,
    with docs/source/scripts/ and creates these directories.

    Parameters
    ----------
    list_fpath : list
        path of scripts to be documented

    Returns
    -------
    list
        path to the directories where documentation (rst files) will be written
    """
    replace_by_path = os.path.join("docs", "source", "scripts")
    list_path_dir_rst = []
    for fpath in list_fpath:
        # extract "/home/siddharth/Autometa/autometa/common/external"
        dir_path = os.path.dirname(os.path.abspath(fpath))
        path_list = dir_path.split("/")
        autometa_index = path_list.index("autometa")
        path_list[autometa_index] = replace_by_path
        # referenced from https://stackoverflow.com/a/14826889/12671809
        # https://docs.python.org/2/tutorial/controlflow.html#unpacking-argument-lists
        path_dir_rst = os.path.join("/", *path_list)
        # path now has "docs/source/scripts/comon/external" instead of autometa
        list_path_dir_rst.append(path_dir_rst)
        if not os.path.exists(path_dir_rst):
            os.makedirs(path_dir_rst)
    return (list_path_dir_rst)


def write_run_argparse_output(list_path_temp, list_path_dir_rst, list_fpath):
    """
    Runs the `--help` command on the agrparse block, copies the output to a rst file
    with proper identation

    Parameters
    ----------
    list_path_temp : list
        temporary files with the argparse block
    list_path_dir_rst : list
        path to the directories where documentation (rst files) will be written
    list_fpath : list
        path of scripts to be documented
    """
    count = 0
    for path_temp in list_path_temp:
        # Extract the filename, eg: kmers.py
        fname = os.path.basename(list_fpath[count])
        # extract the "kmers" of kmers.py
        basename = os.path.splitext(fname)[0]
        fname_rst = basename + ".rst"  # make it kmers.rst
        # another temporary file to store stdout
        fd_temp, path_temp_rst = tempfile.mkstemp()
        cmd = f"python {path_temp} -h"
        # this block writes the argparse output to the respective ".rst" file with proper indentation
        # this is the output that will be under "..code-block::" and will be displayed in html
        with open(path_temp_rst, "w+") as stdout:
            # capture the argparse output
            subprocess.call(cmd, stdout=stdout, shell=True)
            # takes the cursor to beginning to copy and indent stdout.
            stdout.seek(0)
            # The indentation is done because the arparse text needs to be indented
            # under the ..code-block:: section in the final rst file
            # https://docs.python.org/2/library/textwrap.html#textwrap.TextWrapper
            wrapper = textwrap.TextWrapper(
                initial_indent="\t", subsequent_indent="\t", width=80)
            wrapped = ""
            path_script_rst = os.path.join(
                list_path_dir_rst[count], fname_rst)
            count += 1
            # now the path is "docs/source/scripts/comon/kmers.rst"
            with open(path_script_rst, "w") as fh_rst:
                write_text = "\n".join(["="*len(fname), fname, "=" *
                                        len(fname), " ", ".. code-block:: shell", "\n"])
                fh_rst.write(write_text)
                for line in stdout:
                    # copying the indented text
                    wrapped += wrapper.fill(line) + "\n"
                # writing the indented text to the final rst files, this will be dislayed in html
                fh_rst.write(wrapped)
            os.remove(path_temp_rst)


def write_usage_rst():
    """
    Creates the basic outline of usage.rst file
    """
    path_scripts = os.path.join(os.getcwd(), "source", "scripts")
    path_usage = os.path.join(path_scripts, "Usage.rst")
    if not os.path.exists(path_usage):
        with open(path_usage, "w") as fh_main:
            write_text = "\n".join(["="*len("Usage"), "Usage", "="*len("Usage"), " ", ".. toctree::",
                                    "\t :maxdepth: 3", "\t :caption: Table of Contents", " "])
            fh_main.write(write_text)


def design_main_file(list_fpath):
    """
    Makes a list of directories and file name whose toctree will link all the scripts

    Parameters
    ----------
    list_fpath : list
        path of scripts to be documented

    Returns
    -------
    list
        directory name where `_main.rst` files will be and `DirectoryName_main.rst` files
    """

    list_main_file = []
    list_main_dir = []
    for fpath in list_fpath:
        # extract "/home/siddharth/Autometa/autometa/common/external"
        dir_path = os.path.dirname(os.path.abspath(fpath))
        main_dir = os.path.basename(dir_path)  # extraxt "external"
        main_file = main_dir + "_main.rst"  # making a file external_main
        list_main_dir.append(main_dir)
        list_main_file.append(main_file)
    return (list_main_dir, list_main_file)


def design_path_main_file(list_path_dir_rst, list_main_file):
    """
    Makes a list of paths where the `_main.rst` file for each directory will be

    Parameters
    ----------
    list_path_dir_rst : list
        path to the directories where documentation (rst files) will be written
    list_main_file : list
        `DirectoryName_main.rst` files

    Returns
    -------
    list
        path to `DirectoryName_main.rst` files
    """
    count = 0
    list_path_main_file = []
    for path_dir_rst in list_path_dir_rst:
        path_main_file = os.path.join(path_dir_rst, list_main_file[count])
        list_path_main_file.append(path_main_file)
        count += 1
    return (list_path_main_file)


def write_main_rst(list_path_dir_rst, list_main_dir, list_path_main_file):
    """
    Creates a _main file for each directory. These rst files will be called 
    by the toctreee in scripts_main.rst

    Parameters
    ----------
    list_path_dir_rst : list
        path to the directories where documentation (rst files) will be written
    list_main_dir : list
        directory name where `_main.rst` files will be
    list_path_main_file : list
        path to `DirectoryName_main.rst` files
    """
    count = 0
    for path_dir_rst in list_path_dir_rst:
        len_main_dir = len(list_main_dir[count])
        # this will
        if not os.path.exists(list_path_main_file[count]):
            with open(list_path_main_file[count], "w") as fh_main:
                write_text = "\n".join(["="*len_main_dir, list_main_dir[count], "="*len_main_dir, " ", ".. toctree::",
                                        "\t :maxdepth: 2", "\t :caption: Table of Contents", " "])
                fh_main.write(write_text)
        count += 1


def link_sub_dir(list_fpath, list_path_main_file):
    """
    Links sub-directories with the toc tree of parent directory.
    Eg. Linking external with the toctree of common

    Parameters
    ----------
    list_fpath : list
        path of scripts to be documented
    list_path_main_file : list
        path to `DirectoryName_main.rst` files
    """
    count = 0
    path_recorded = ""
    for fpath in list_fpath:
        dir_path = os.path.dirname(os.path.abspath(fpath))
        for root, dirnames, filename in os.walk(dir_path):
            # This is done to prevent writing directories below  the sub-directories
            # i.e. we don't need to link any directories (if they are in future) that are below external with the common toctree
            # Done to take into account any future upgrades
            if path_recorded == list_path_main_file[count]:
                continue
            # this counts the number of sub-directories each directory has
            # true only for `common`
            if len(dirnames) > 0:
                for dirname in dirnames:
                    dir_main = dirname + "_main"  # external_main
                    path_recorded = list_path_main_file[count]
                    with open(list_path_main_file[count], "a+") as fh_main:
                        fh_main.write(f"\n\t {dirname}/{dir_main}")
                        # writes external/external_main in common_main
    count += 1


def write_script_name(list_path_main_file, list_fpath):
    """
    Opens each DirectoryName_main.rst and add the name of the scripts in that directory
    # Eg. in external.rst, we add samtool, bedtools, prodigal, etc

    Parameters
    ----------
    list_path_main_file : list
        path to `DirectoryName_main.rst` files
    list_fpath : list
        path of scripts to be documented
    """
    count = 0
    for path_main_file in list_path_main_file:
        fname = os.path.basename(list_fpath[count])
        # extract the "kmers" of kmers.py
        basename = os.path.splitext(fname)[0]
        with open(path_main_file, "a") as fh_main:
            fh_main.write(f"\n\t {basename}")
        count += 1


def link_dir():
    """
    Writes the directories in usage.rst toctree
    """
    # This can also be done using next() -> https://stackoverflow.com/a/142535/12671809 , but reduces readability. What to do ??
    path_scripts = os.path.join(os.getcwd(), "source", "scripts")
    path_usage = os.path.join(path_scripts, "Usage.rst")
    for dirname in os.listdir(path_scripts):
        path_dirname = os.path.join(path_scripts, dirname)
        if os.path.isdir(path_dirname):
            with open(path_usage, "a") as fh_usage:
                dirname_rst = dirname + "_main"
                fh_usage.write(f"\n\t {dirname}/{dirname_rst}")


def remove_existing_docs():
    """
    Removes the exists `scripts` directory
    """
    path_scripts = os.path.join(os.getcwd(), "source", "scripts")
    if os.path.exists(path_scripts):
        shutil.rmtree(path_scripts)


def remove_empty_dir():
    """
    Removes any empty directory that may have been added by mistake
    """
    for root, dirnames, filenames in os.walk(os.path.join(os.path.dirname(os.getcwd()), "autometa")):
        if os.listdir(root) == []:
            os.rmdir(root)


def remove_temp_files(list_path_temp):
    """
    Removes the temporary files to which the argparse block was copied

    Parameters
    ----------
    list_path_temp : list
        temporary files with the argparse block
    """
    for i in list_path_temp:
        os.remove(i)


remove_existing_docs()
remove_empty_dir()

# get the list of files which needs to be documented
list_fpath = parse_files()
# copy argparse block to temp files
list_path_temp = get_argparse_block(list_fpath)
# list of directory paths where rst files will be written and create those directories
list_path_dir_rst = path_dir_docs_scripts(list_fpath)
# run the argparse block and copy the output to respective rst files
write_run_argparse_output(list_path_temp, list_path_dir_rst, list_fpath)
# Makes a list of directories and file name whose toctree will link all the scripts
list_main_dir, list_main_file = design_main_file(list_fpath)
# Makes a list of paths where the `_main.rst` file for each directory will be
list_path_main_file = design_path_main_file(list_path_dir_rst, list_main_file)
# Creates a _main file for each directory
write_main_rst(list_path_dir_rst, list_main_dir, list_path_main_file)
# Links sub-directories with the toc tree of parent directory.
link_sub_dir(list_fpath, list_path_main_file)
# Opens each DirectoryName_main.rst and add the name of the scripts in that directory
write_script_name(list_path_main_file, list_fpath)
# create the top of usage.rst file
write_usage_rst()
# Writes the directories in usage.rst toctree
link_dir()

remove_temp_files(list_path_temp)
