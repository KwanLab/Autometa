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
import os
import re
import shutil
import subprocess
import tempfile
import textwrap


def parse_files():
    """
    Makes a list of all the scripts that needs to be documented

    Returns
    -------
    list
        file paths that needs to be documented, 
        e.g. ['</path/to/autometa/script.py>', ...]
    """
    exclude_dirs = ["validation"]
    fpaths = []
    for fpath in glob.glob(os.path.join(os.path.dirname(os.path.dirname(__file__)), "autometa", "**", "*.py"), recursive=True):
        if "__" in fpath:  # do not include __init__ and __main__.py files
            continue
        if os.path.basename(os.path.dirname(fpath)) in exclude_dirs:
            continue
        fpaths.append(fpath)
    return(fpaths)


def get_argparse_block(fpaths):
    """
    Copies the argparse block of each script in a temporary file.

    Parameters
    ----------
    fpaths : list
        file paths that needs to be documented,
        e.g. ['</path/to/autometa/script.py>', ...]

    Returns
    -------
    list
        temporary files with the argparse block
        e.g. ['</path/to/tmp/file>', ...]
    """
    tmp_fpaths = []
    for fpath in fpaths:
        writing = False
        outlines = ""
        with open(fpath) as fh:
            for line in fh:
                line = line.strip()
                if line == "import argparse":
                    writing = True
                    imports = ["import os",
                               "import multiprocessing as mp", "\n"]
                    import_lines = "\n".join(imports)
                    outlines += import_lines
                # Add usage  = ScriptName.py in argparse block
                if "argparse.ArgumentParser" in line:
                    fname = os.path.basename(fpath)
                    # Uses the starting `(` of Argument.parser, and replaces it with ( usage=
                    usage = f"( usage = \"{fname}\", "
                    line = line.replace("(", usage)
                if writing:  # Convert default constants to `str`
                    # eg. config.DEFAULT_FPATH to "config.DEFAULT_FPATH"
                    x = re.search(r"[a-z]*\.?[A-Z]+_[A-Z]+[_A-Z]*", line)
                    if x:
                        line = line.replace(x.group(), f'"{x.group()}"')
                if writing:
                    outlines += f"{line}\n"
                    if line == "args = parser.parse_args()":
                        writing = False
                __, tmp_fpath = tempfile.mkstemp()
                with open(tmp_fpath, 'w') as outfh:
                    outfh.write(outlines)
        tmp_fpaths.append(tmp_fpath)
    return tmp_fpaths


def path_dir_docs_scripts(fpaths):
    """
    Replaces the autometa directory in the path of scripts that needs to be documented,
    with <docs/source/scripts/> and creates these directories.

    Parameters
    ----------
    fpaths : list
        file paths that needs to be documented, 
        e.g. ['</path/to/autometa/script.py>', ...]

    Returns
    -------
    list
        path to the directories where directories where script.rst files will be written.
        e.g. ['</path/to/docs/source/scripts/directory>', ...]
    """
    rst_scripts_dir = os.path.join("docs", "source", "scripts")
    rst_dirpaths = []
    for fpath in fpaths:
        # change this path "~/Autometa/autometa/common/external"
        script_dirpath = os.path.dirname(os.path.abspath(fpath))
        path_list = script_dirpath.split(os.path.sep)
        autometa_index = path_list.index("autometa")
        path_list[autometa_index] = rst_scripts_dir
        # referenced from https://stackoverflow.com/a/14826889/12671809
        # https://docs.python.org/2/tutorial/controlflow.html#unpacking-argument-lists
        rst_dirpath = os.path.join(os.path.sep, *path_list)
        # path now has "~/Autometa/docs/source/scripts/comon/external" instead of autometa
        if not os.path.exists(rst_dirpath):
            os.makedirs(rst_dirpath)
        rst_dirpaths.append(rst_dirpath)
    return (rst_dirpaths)


def write_run_argparse_output(tmp_fpaths, rst_dirpaths, fpaths):
    """
    Runs the `--help` command on the agrparse block, copies the output to a rst file
    with proper identation

    Parameters
    ----------
    tmp_fpaths : list
        temporary files with the argparse block
        e.g. ['</path/to/tmp/file>', ...]
    rst_dirpaths : list
        path to the directories where script.rst files will be written.
        e.g. ['</path/to/docs/source/scripts/directory>', ...]
    fpaths : list
        file paths that needs to be documented, 
        e.g. ['</path/to/autometa/script.py>', ...]
    """
    count = 0
    for path_temp in tmp_fpaths:
        # Extract the filename, eg: kmers.py
        fname = os.path.basename(fpaths[count])
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
                rst_dirpaths[count], fname_rst)
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


def design_main_file(fpaths):
    """
    Makes a list of directories and file name whose toctree will link all the scripts

    Parameters
    ----------
    fpaths : list
        file paths that needs to be documented, 
        e.g. ['</path/to/autometa/script.py>', ...]

    Returns
    -------
    list
        directory name where `_main.rst` files will be and `DirectoryName_main.rst` files
    """

    rst_main_files = []
    rst_main_dirs = []
    for fpath in fpaths:
        # extract "/home/siddharth/Autometa/autometa/common/external"
        dir_path = os.path.dirname(os.path.abspath(fpath))
        rst_main_dir = os.path.basename(dir_path)  # extraxt "external"
        rst_main_file = rst_main_dir + "_main.rst"  # making a file external_main
        rst_main_dirs.append(rst_main_dir)
        rst_main_files.append(rst_main_file)
    return (rst_main_dirs, rst_main_files)


def design_path_main_file(rst_dirpaths, rst_main_files):
    """
    Makes a list of paths where the `_main.rst` file for each directory will be

    Parameters
    ----------
    rst_dirpaths : list
        path to the directories where directories where script.rst files will be written.
        e.g. ['</path/to/docs/source/scripts/directory>', ...]
    rst_main_files : list
        `DirectoryName_main.rst` files

    Returns
    -------
    list
        path to `DirectoryName_main.rst` files
    """
    count = 0
    path_main_files = []
    for rst_dirpath in rst_dirpaths:
        path_main_file = os.path.join(rst_dirpath, rst_main_files[count])
        path_main_files.append(path_main_file)
        count += 1
    return (path_main_files)


def write_main_rst(rst_dirpaths, rst_main_dirs, path_main_files):
    """
    Creates a _main file for each directory. These rst files will be called 
    by the toctreee in scripts_main.rst

    Parameters
    ----------
    rst_dirpaths : list
        path to the directories where directories where script.rst files will be written.
        e.g. ['</path/to/docs/source/scripts/directory>', ...]
    rst_main_dirs : list
        directory name where `_main.rst` files will be
    path_main_files : list
        path to `DirectoryName_main.rst` files
    """
    count = 0
    for rst_dirpath in rst_dirpaths:
        len_main_dir = len(rst_main_dirs[count])
        if not os.path.exists(path_main_files[count]):
            with open(path_main_files[count], "w") as fh_main:
                write_text = "\n".join(["="*len_main_dir, rst_main_dirs[count], "="*len_main_dir, " ", ".. toctree::",
                                        "\t :maxdepth: 2", "\t :caption: Table of Contents", " "])
                fh_main.write(write_text)
        count += 1


def link_sub_dir(fpaths, path_main_files):
    """
    Links sub-directories with the toc tree of parent directory.
    Eg. Linking external with the toctree of common

    Parameters
    ----------
    fpaths : list
        file paths that needs to be documented, 
        e.g. ['</path/to/autometa/script.py>', ...]
    path_main_files : list
        path to `DirectoryName_main.rst` files
    """
    count = 0
    path_recorded = ""
    for fpath in fpaths:
        dir_path = os.path.dirname(os.path.abspath(fpath))
        # print(dir_path)
        for root, dirnames, filename in os.walk(dir_path):
            # This is done to prevent writing directories below  the sub-directories
            # i.e. we don't need to link any directories (if they are in future) that are below external with the common toctree
            # Done to take into account any future upgrades
            if path_recorded == path_main_files[count]:
                continue
            # this counts the number of sub-directories each directory has
            # true only for `common`
            if len(dirnames) > 0:
                for dirname in dirnames:
                    if "__" in dirname:
                        continue
                    dir_main = dirname + "_main"  # external_main
                    path_recorded = path_main_files[count]
                    with open(path_main_files[count], "a+") as fh_main:
                        fh_main.write(f"\n\t {dirname}/{dir_main}")
                        # writes external/external_main in common_main
        count += 1


def write_script_name(path_main_files, fpaths):
    """
    Opens each DirectoryName_main.rst and add the name of the scripts in that directory
    # Eg. in external.rst, we add samtool, bedtools, prodigal, etc

    Parameters
    ----------
    path_main_files : list
        path to `DirectoryName_main.rst` files
    fpaths : list
        file paths that needs to be documented, 
        e.g. ['</path/to/autometa/script.py>', ...]
    """
    count = 0
    for path_main_file in path_main_files:
        fname = os.path.basename(fpaths[count])
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


def remove_temp_files(tmp_fpaths):
    """
    Removes the temporary files to which the argparse block was copied

    Parameters
    ----------
    tmp_fpaths : list
        temporary files with the argparse block
        e.g. ['</path/to/tmp/file>', ...]
    """
    for i in tmp_fpaths:
        os.remove(i)


remove_existing_docs()
remove_empty_dir()

# get the list of files which needs to be documented
fpaths = parse_files()
# copy argparse block to temp files
tmp_fpaths = get_argparse_block(fpaths)
# list of directory paths where rst files will be written and create those directories
rst_dirpaths = path_dir_docs_scripts(fpaths)
# run the argparse block and copy the output to respective rst files
write_run_argparse_output(tmp_fpaths, rst_dirpaths, fpaths)
# Makes a list of directories and file name whose toctree will link all the scripts
rst_main_dirs, rst_main_files = design_main_file(fpaths)
# Makes a list of paths where the `_main.rst` file for each directory will be
path_main_files = design_path_main_file(rst_dirpaths, rst_main_files)
# Creates a _main file for each directory
write_main_rst(rst_dirpaths, rst_main_dirs, path_main_files)
# Links sub-directories with the toc tree of parent directory.
link_sub_dir(fpaths, path_main_files)
# Opens each DirectoryName_main.rst and add the name of the scripts in that directory
write_script_name(path_main_files, fpaths)
# create the top of usage.rst file
write_usage_rst()
# Writes the directories in usage.rst toctree
link_dir()

remove_temp_files(tmp_fpaths)
