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


def get_argparse_block(fpath):
    """
    Copies the argparse block of each script in a string

    Regex is used to capture the constant string eg. config.DEFAULT_FPATH then convert it to a string eg. "config.DEFAULT_FPATH"
    Regex functioning: requied 'default=', lower case letter (optional) followed by a '.' (zero or one) upper case (atleast one)
    the optional '_' (to match something like DEFAULT), upper case (atleast one), then underscore ('_') followed by upper case, both os them optional
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
            match = re.search(
                r"default=([a-z]*\.?[A-Z]+_?[A-Z]+[_A-Z]*)", line)
            if match:
                cap_const = match.group(1)
                line = line.replace(
                    cap_const, f'"{cap_const}"')
            outlines += f"{line}\n"
    return outlines


def get_uasge(argparse_lines, fpath):
    """
    Write the argparse block  to a temporary file and run the `--help` command on it,
    followed by proper indentation as per rst syntax

    Parameters
    ----------
    argparse_lines : string
        lines returned from `get_argparse_block` i.e. argparse block between (and includes) import argparse and
        args = parser.parse_args() to be written to a temporary file
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'

    Returns
    -------
    wrapped_lines : string
        indented arparse output after running the `--help` command
    """
    __, tmp_fpath = tempfile.mkstemp()
    with open(tmp_fpath, 'w') as outfh:
        outfh.write(argparse_lines)
    cmd = f"python {tmp_fpath} -h"
    # capture the argparse output
    proc = subprocess.run(cmd, stdout=subprocess.PIPE,
                          stderr=subprocess.DEVNULL, shell=True, text=True)
    # checks if the command was run properly else raises error
    try:
        proc.check_returncode()
    except:
        print(f"Error while parsing argparse block of {fpath}")
    wrapper = textwrap.TextWrapper(
        initial_indent="\t", subsequent_indent="\t", width=80)
    wrapped_lines = ""
    # splitlines is required because wrapped_lines.stdout is just as big line with \n in between.
    # Splitlines will "split" as each \n, without this var `line` would go through each character and not each line
    for line in proc.stdout.splitlines():
        # writing the indented text to the final rst files, this will be dislayed in html
        wrapped_lines += wrapper.fill(line) + "\n"
    os.remove(tmp_fpath)
    return (wrapped_lines)


def make_rst_dir(fpath):
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


def write_usage(fpath, rst_dirpath, wrapped_lines):
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
    wrapped_lines : string
        indented arparse output after running the `--help` command
    """
    # Extract the filename, eg: kmers.py
    fname = os.path.basename(fpath)
    # extract the "kmers" of kmers.py
    basename = os.path.splitext(fname)[0]
    fname_rst = basename + ".rst"  # make it kmers.rst
    path_script_rst = os.path.join(
        rst_dirpath, fname_rst)
    frmt_header = "="*len(fname)
    text = f"{frmt_header}\n{fname}\n{frmt_header}\n\n.. code-block:: shell\n"
    with open(path_script_rst, "w") as fh_rst:
        fh_rst.write(f"{text}\n{wrapped_lines}")


def write_main_file(fpath, rst_dirpath):
    """
    Writes the header and basic information to all of the `main.rst` files

    Parameters
    ----------
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'
    rst_dirpath : string
        path to the directory where script.rst files will be written,
        e.g. '</path/to/docs/source/scripts/directory>'
    """
    # extract "/home/siddharth/Autometa/autometa/common/external"
    dir_path = os.path.dirname(os.path.abspath(fpath))
    main_rst_dir = os.path.basename(dir_path)  # extraxt "external"
    main_rst_file = main_rst_dir + "_main.rst"  # making a file external_main
    main_rst_fpath = os.path.join(rst_dirpath, main_rst_file)
    if not os.path.exists(main_rst_fpath):
        frmt_header = "="*len(main_rst_dir)
        outlines = f"{frmt_header}\n{main_rst_dir}\n{frmt_header}\n\n.. toctree::\n"
        outlines += f"\t :maxdepth: 2\n\t :caption: Table of Contents\n\n"
        with open(main_rst_fpath, "w") as fh:
            fh.write(outlines)


def link_main(rst_scripts_dirpath):
    """
    Links the parent toc tree with the sub-directories and the script that
    directory has

    Parameters
    ----------
    rst_scripts_dirpath : string
        path to the directory where script.rst files will be written,
        e.g. '</path/to/docs/source/scripts/directory>'
    """
    for root, dirnames, filenames in os.walk(rst_scripts_dirpath):
        if root == rst_scripts_dirpath:
            continue
        # links sub-directories
        lines = ""
        if dirnames:
            for dirname in dirnames:
                dir_main = dirname + "_main"  # external_main
                lines += f"\n\t {dirname}/{dir_main}"
        for index, filename in enumerate(filenames):
            # capture the index of `_main.rst` to open later
            if "main" in filename:
                break
        # get the path of the `main.rst` file
        main_rst_fpath = os.path.join(root, filenames[index])
        # write the name of each script in the `main.rst` file
        for fname in filenames:
            # skips the `main.rst` as it is not an Autometa script, and
            # it is the file we are writing to
            if "main" in fname:
                continue
            fname = os.path.splitext(fname)[0]
            lines += f"\n\t {fname}"
        with open(main_rst_fpath, "a") as fh:
            fh.write(lines)


def write_usage_file(rst_scripts_dirpath):
    """
    Creates usage.rst and writes the directories in its toctree

    Parameters
    ----------
    rst_scripts_dirpath : string
        path to the directory where script.rst files will be written,
        e.g. '</path/to/docs/source/scripts/directory>'
    """
    usage_fpath = os.path.join(rst_scripts_dirpath, "usage.rst")
    frmt_header = "="*len("usage")
    text = f"{frmt_header}\nUsage\n{frmt_header}\n\n.. toctree::\n"
    text += f"\t :maxdepth: 3\n\t :caption: Table of Contents\n\n"
    if not os.path.exists(usage_fpath):
        with open(usage_fpath, "w") as fh:
            fh.write(text)
            for dirname in os.listdir(rst_scripts_dirpath):
                dirpath = os.path.join(rst_scripts_dirpath, dirname)
                if os.path.isdir(dirpath):
                    dirname_rst = dirname + "_main"
                    fh.write(f"\n\t {dirname}/{dirname_rst}")


def remove_existing_docs(rst_scripts_dirpath):
    """
    Removes the existing documentatio in `scripts` directory

    Parameters
    ----------
    rst_scripts_dirpath : string
        path to the directory where script.rst files will be written,
        e.g. '</path/to/docs/source/scripts/directory>'
    """
    if os.path.exists(rst_scripts_dirpath):
        shutil.rmtree(rst_scripts_dirpath)


def remove_empty_dir():
    """
    Removes any empty directory that may have been added by mistake
    """
    pkg_dirpath = os.path.join(os.path.dirname(
        os.path.dirname(__file__)), "autometa")
    for root, _, __ in os.walk(pkg_dirpath):
        # If the list is empty, then this will evaluate to false, and if it contains elements, it will evaluate to true
        if not os.listdir(root):
            os.rmdir(root)


def main():
    rst_scripts_dirpath = os.path.join(
        os.path.dirname(__file__), "source", "scripts")
    remove_existing_docs(rst_scripts_dirpath)
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
        wrapped_lines = get_uasge(argparse_lines, fpath)
        rst_dirpath = make_rst_dir(fpath)
        write_usage(fpath, rst_dirpath, wrapped_lines)
        write_main_file(fpath, rst_dirpath)
    link_main(rst_scripts_dirpath)
    write_usage_file(rst_scripts_dirpath)


if __name__ == '__main__':
    main()
