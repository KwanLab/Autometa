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


def get_argparse_block(fpath):
    """
    Copies the argparse block of each script in a string

    Regex is used to capture the constant string eg. config.DEFAULT_FPATH then convert it to a string eg. "config.DEFAULT_FPATH"
    Regex functioning: lower case letter (optional) followed by a '.' (zero or one) upper case (atleast one)
    the required '_', upper case (atleast one), then underscore ('_') followed by upper case, both os them optional
    captures, config.DEFAULT_FPATH, MARKERS_DIR, wq.WORK_QUEUE_DEFAULT_PORT and DEFAULT_FPATH

    Parameters
    ----------
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'

    Returns
    -------
    string
        argparse block between (and includes) import argparse and args = parser.parse_args() to be written to a temporary file
    """
    writing = False
    outlines = "import os\nimport multiprocessing as mp\n"
    with open(fpath) as fh:
        for line in fh:
            line = line.strip()
            if line == "import argparse":
                writing = True
            if line == "args = parser.parse_args()":
                outlines += f"{line}\n"
                writing = False
            if not writing:
                continue
            # Notice below we can assume we are writing and no longer need additional if statements
            # Add usage  = script.py in argparse block
            if "argparse.ArgumentParser" in line:
                fname = os.path.basename(fpath)
                # Uses the starting `(` of Argument.parser, and replaces it with ( usage=
                # this adds the usage = script.py in Argparse automatically
                usage = f"(usage='{fname}', "
                line = line.replace("(", usage)
            # Regex is used to capture the constant string eg. config.DEFAULT_FPATH then convert it to a string eg. "config.DEFAULT_FPATH"
            capture_constant = re.search(
                r"[a-z]*\.?[A-Z]+_[A-Z]+[_A-Z]*", line)
            if capture_constant:
                line = line.replace(
                    capture_constant.group(), f'"{capture_constant.group()}"')
            outlines += f"{line}\n"
    return outlines


def path_dir_docs_scripts(fpath):
    """
    Replaces the autometa directory in the path of scripts that needs to be documented,
    with <docs/source/scripts/> and creates these directories

    Parameters
    ----------
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'

    Returns
    -------
    string
        path to the directory where script.rst files will be written,
        e.g. '</path/to/docs/source/scripts/directory>'
    """
    rst_scripts_dir = os.path.join("docs", "source", "scripts")
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
    return (rst_dirpath)


def get_argparse_out(argparse_lines):
    """
    Write the argparse block to a temporary file and run the `--help` command on it
    Indents the output using `sed` to be inaccordance with rst syntax

    Parameters
    ----------
    argparse_lines : string
        argparse block between (and includes) import argparse and 
        args = parser.parse_args() to be written to a temporary file

    Returns
    -------
    tmp_fpath : string
        path to temporary file that has the argparse block
        e.g. '</path/to/tmp/file>'
    proc : string
        indented arparse output after running the `--help` command    
    """
    __, tmp_fpath = tempfile.mkstemp()
    with open(tmp_fpath, 'w') as outfh:
        outfh.write(argparse_lines)
    cmd = f"python {tmp_fpath} -h | sed 's/^/   /'"
    # capture the argparse output
    proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, shell=True, text=True)
    return (tmp_fpath, proc)


def write_argparse_out(fpath, rst_dirpath, proc):
    """
    Creates `rst` file for each script, and writes the header along with 
    the respective argparse output generated after running `--help`

    Parameters
    ----------
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'
    rst_dirpath : string
        path to the directory where script.rst files will be written,
        e.g. '</path/to/docs/source/scripts/directory>'
    proc : string
        indented arparse output after running the `--help` command    
    """
    # Extract the filename, eg: kmers.py
    fname = os.path.basename(fpath)
    # extract the "kmers" of kmers.py
    basename = os.path.splitext(fname)[0]
    fname_rst = basename + ".rst"  # make it kmers.rst
    path_script_rst = os.path.join(
        rst_dirpath, fname_rst)
    with open(path_script_rst, "w") as fh_rst:
        frmt_header = "="*len(fname)
        text = f"{frmt_header}\n{fname}\n{frmt_header}\n\n.. code-block:: shell\n"
        fh_rst.write(f"{text}\n{proc.stdout}")


def design_main_file(fpath, rst_dirpath):
    """
    Creates a _main file for each directory, and link sub-directories with the toc tree of parent directory.
    Eg. Linking external with the toctree of common
    These rst files will be called by the toctreee in scripts_main.rst

    Also prevents writing directories below  the sub-directories with the toc tree of parent directory,
    ie. we don't need to link any directories (if they are in future) that are below external with the common toctree.
    The way this `if` statements works is: os.walk will keep iterating through a directory recursively (going into sub-sub-directories), we just need its first output
    of the directories within the root, after that it goes through the sub-sub-directories and lists their directories which we don't want,
    meanwhile main_file_fpath is constant, thus we can skip recording sub-subdirectories if we tell it that, this path has alreasy been recorded.
    Remember main_file_fpath is remaining constant for one dir_path


    Parameters
    ----------
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'
    rst_dirpath : string
        path to the directory where script.rst files will be written,
        e.g. '</path/to/docs/source/scripts/directory>'

    Returns
    -------
    string
        path to the `_main.rst` file corresponding to each directory
        e.g. '</path/to/docs/source/scripts/directory/directory_main.rst>' 
    """
    # extract "/home/siddharth/Autometa/autometa/common/external"
    dir_path = os.path.dirname(os.path.abspath(fpath))
    rst_main_dir = os.path.basename(dir_path)  # extraxt "external"
    rst_main_file = rst_main_dir + "_main.rst"  # making a file external_main
    main_file_fpath = os.path.join(rst_dirpath, rst_main_file)
    outlines = ""
    if os.path.exists(main_file_fpath):
        return main_file_fpath
    with open(main_file_fpath, "w") as fh:
        frmt_header = "="*len(rst_main_dir)
        outlines += f"{frmt_header}\n{rst_main_dir}\n{frmt_header}\n\n.. toctree::\n"
        outlines += f"\t :maxdepth: 2\n\t :caption: Table of Contents\n\n"
        path_recorded = ""
        for __, dirnames, _ in os.walk(dir_path):
            if path_recorded == main_file_fpath:
                continue
            for dirname in dirnames:
                if "__" in dirname:  # to prevent __pychache__ folders from being considered
                    continue
                dir_main = dirname + "_main"  # external_main
                path_recorded = main_file_fpath
                lines = f"\n\t {dirname}/{dir_main}"
                outlines += lines
        fh.write(outlines)
    return (main_file_fpath)


def write_scrit_names(main_file_fpath, fpath):
    """
    Opens each DirectoryName_main.rst and add the name of the scripts in that directory

    Parameters
    ----------
    main_file_fpath : string
        path to the `_main.rst` file corresponding to each directory
        e.g. '</path/to/docs/source/scripts/directory/directory_main.rst>' 
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'
    """
    with open(main_file_fpath, "a") as fh:
        basename = os.path.basename(fpath)
        # extract the "kmers" of kmers.py
        fname = os.path.splitext(basename)[0]
        fh.write(f"\n\t {fname}")


def link_write_usage():
    """
    Creates usage.rst and writes the directories in its toctree
    """
    scripts_dirpath = os.path.join(
        os.path.dirname(__file__), "source", "scripts")
    usage_fpath = os.path.join(scripts_dirpath, "usage.rst")
    if not os.path.exists(usage_fpath):
        with open(usage_fpath, "w") as fh:
            frmt_header = "="*len("usage")
            text = f"{frmt_header}\nUsage\n{frmt_header}\n\n.. toctree::\n"
            text += f"\t :maxdepth: 3\n\t :caption: Table of Contents\n\n"
            fh.write(text)
            for dirname in os.listdir(scripts_dirpath):
                dirpath = os.path.join(scripts_dirpath, dirname)
                if os.path.isdir(dirpath):
                    dirname_rst = dirname + "_main"
                    fh.write(f"\n\t {dirname}/{dirname_rst}")


def remove_existing_docs():
    """
    Removes the existing documentatio in `scripts` directory
    """
    scripts_dirpath = os.path.join(
        os.path.dirname(__file__), "source", "scripts")
    if os.path.exists(scripts_dirpath):
        shutil.rmtree(scripts_dirpath)


def remove_empty_dir():
    """
    Removes any empty directory that may have been added by mistake
    """
    pkg_dirpath = os.path.join(os.path.dirname(
        os.path.dirname(__file__)), "autometa")
    for root, dirnames, filenames in os.walk(pkg_dirpath):
        # If the list is empty, then this will evaluate to false, and if it contains elements, it will evaluate to true
        if not os.listdir(root):
            os.rmdir(root)


def main():
    remove_existing_docs()
    remove_empty_dir()
    pkg_scripts_fpath = os.path.join(os.path.dirname(
        os.path.dirname(__file__)), "autometa", "**", "*.py")
    for fpath in glob.glob(pkg_scripts_fpath, recursive=True):
        exclude_dirs = ["validation"]
        if "__" in fpath:  # do not include __init__ and __main__.py files
            continue
        if os.path.basename(os.path.dirname(fpath)) in exclude_dirs:
            continue
        argparse_lines = get_argparse_block(fpath)
        rst_dirpath = path_dir_docs_scripts(fpath)
        tmp_fpath, proc = get_argparse_out(argparse_lines)
        write_argparse_out(fpath, rst_dirpath, proc)
        main_file_fpath = design_main_file(fpath, rst_dirpath)
        write_scrit_names(main_file_fpath, fpath)
    link_write_usage()
    os.remove(tmp_fpath)


if __name__ == '__main__':
    main()
