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
    Copies the argparse block of each script in a temporary file

    Parameters
    ----------
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'

    Returns
    -------
    string
        path to temporary file that has the argparse block
        e.g. '</path/to/tmp/file>'
    """
    writing = False
    outlines = ""
    with open(fpath) as fh:
        for line in fh:
            line = line.strip()
            if line == "import argparse":
                writing = True
                imports = ["import os",
                           "import multiprocessing as mp", ""]
                # the blank "" is needed to give a new line else error
                import_lines = "\n".join(imports)
                outlines += import_lines
            # Add usage  = ScriptName.py in argparse block
            if "argparse.ArgumentParser" in line:
                fname = os.path.basename(fpath)
                # Uses the starting `(` of Argument.parser, and replaces it with ( usage=
                usage = f"( usage = '{fname}', "
                line = line.replace("(", usage)
            if writing:
                # regex is used to capture the constant string eg. config.DEFAULT_FPATH
                # then convert it to "config.DEFAULT_FPATH"
                # Regex functioning: lower case letter (optional) followed by a '.' (zero or one) upper case (atleast one)
                # the required '_', upper case (atleast one), then underscore ('_') followed by upper case, both os them optional
                #captures, config.DEFAULT_FPATH, MARKERS_DIR, wq.WORK_QUEUE_DEFAULT_PORT and DEFAULT_FPATH
                # add the a comment there and change the variable 'cap_str'
                cap_str = re.search(r"[a-z]*\.?[A-Z]+_[A-Z]+[_A-Z]*", line)
                if cap_str:
                    line = line.replace(
                        cap_str.group(), f'"{cap_str.group()}"')
            if writing:
                outlines += f"{line}\n"
            if line == "args = parser.parse_args()":
                writing = False
    __, tmp_fpath = tempfile.mkstemp()
    with open(tmp_fpath, 'w') as outfh:
        outfh.write(outlines)
    return tmp_fpath


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


def write_run_argparse_output(tmp_fpath, rst_dirpath, fpath):
    """
    Runs the `--help` command on the agrparse block, copies the output to a rst file
    with proper identation

    Parameters
    ----------
    tmp_fpath : string
        path to temporary file that has the argparse block
        e.g. '</path/to/tmp/file>'
    rst_dirpath : string
        path to the directory where script.rst files will be written,
        e.g. '</path/to/docs/source/scripts/directory>'
    fpath : string
        file path that needs to be documented,
        e.g. '</path/to/autometa/script.py>'
    """
    # Extract the filename, eg: kmers.py
    fname = os.path.basename(fpath)
    # extract the "kmers" of kmers.py
    basename = os.path.splitext(fname)[0]
    fname_rst = basename + ".rst"  # make it kmers.rst
    # another temporary file to store stdout
    __, path_tmp_rst = tempfile.mkstemp()
    cmd = f"python {tmp_fpath} -h"
    # this block writes the argparse output to the respective ".rst" file with proper indentation
    # this is the output that will be under "..code-block::" and will be displayed in html
    with open(path_tmp_rst, "w+") as stdout:
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
            rst_dirpath, fname_rst)
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
        os.remove(path_tmp_rst)


def design_main_file(fpath, rst_dirpath):
    """
    Creates a _main file for each directory, and link sub-directories with the toc tree of parent directory.
    Eg. Linking external with the toctree of common
    These rst files will be called by the toctreee in scripts_main.rst

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
    rst_main_dir = os.path.basename(dir_path)  # extraxt "external"
    rst_main_file = rst_main_dir + "_main.rst"  # making a file external_main
    path_main_file = os.path.join(rst_dirpath, rst_main_file)
    outlines = ""
    if not os.path.exists(path_main_file):
        with open(path_main_file, "w") as fh:
            lines = "\n".join(["="*len(rst_main_dir), rst_main_dir, "="*len(rst_main_dir), " ", ".. toctree::",
                               "\t :maxdepth: 2", "\t :caption: Table of Contents", " "])
            outlines += lines
            path_recorded = ""
            for __, dirnames, filename in os.walk(dir_path):
                # This is done to prevent writing directories below  the sub-directories
                # i.e. we don't need to link any directories (if they are in future) that are below external with the common toctree
                # the way this `if` statements works is: os.walk will keep iterating through a directory recursively (going into sub-sub-directories), we just need its first output
                # of the directories within the root, after that it goes through the sub-sub-directories and lists their directories which we don't want,
                # meanwhile path_main_file is constant, thus we can skip recording sub-subdirectories if we tell it that, this path has alreasy been recorded.
                # remember path_main_file is remaining constant for one dir_path
                if path_recorded == path_main_file:
                    continue
                # this counts the number of sub-directories each directory has
                # true only for `common` (as of April 24, 2020)
                if len(dirnames) > 0:
                    for dirname in dirnames:
                        if "__" in dirname:  # to prevent __pychache__ folders from being considered
                            continue
                        dir_main = dirname + "_main"  # external_main
                        path_recorded = path_main_file
                        lines = f"\n\t {dirname}/{dir_main}"
                        outlines += lines
            fh.write(outlines)
    # Opens each DirectoryName_main.rst and add the name of the scripts in that directory
    with open(path_main_file, "a") as fh:
        fname = os.path.basename(fpath)
        # extract the "kmers" of kmers.py
        basename = os.path.splitext(fname)[0]
        fh.write(f"\n\t {basename}")


def link_write_usage():
    """
    Creates usage.rst and writes the directories in its toctree
    """
    if not os.path.exists(path_usage):
        with open(path_usage, "w") as fh:
            write_text = "\n".join(["="*len("usage"), "Usage", "="*len("usage"), " ", ".. toctree::",
                                    "\t :maxdepth: 3", "\t :caption: Table of Contents", " "])
            fh.write(write_text)
            # This can also be done using next() -> https://stackoverflow.com/a/142535/12671809 ,
            # but reduces readability. What to do ??
            for dirname in os.listdir(path_scripts):
                path_dirname = os.path.join(path_scripts, dirname)
                if os.path.isdir(path_dirname):
                    dirname_rst = dirname + "_main"
                    fh.write(f"\n\t {dirname}/{dirname_rst}")


def remove_existing_docs():
    """
    Removes the exists documentatio in `scripts` directory
    """
    if os.path.exists(path_scripts):
        shutil.rmtree(path_scripts)


def remove_empty_dir():
    """
    Removes any empty directory that may have been added by mistake
    """
    for root, dirnames, filenames in os.walk(os.path.join(os.path.dirname(os.path.dirname(__file__)), "autometa")):
        if os.listdir(root) == []:
            os.rmdir(root)


def remove_temp_files(tmp_fpath):
    """
    Removes the temporary files to which the argparse block was copied

    Parameters
    ----------
    tmp_fpath : string
        path to temporary file that has the argparse block
        e.g. '</path/to/tmp/file>'
    """
    os.remove(tmp_fpath)


path_scripts = os.path.join(os.path.dirname(__file__), "source", "scripts")
path_usage = os.path.join(path_scripts, "usage.rst")

remove_existing_docs()
remove_empty_dir()

for fpath in glob.glob(os.path.join(os.path.dirname(os.path.dirname(__file__)), "autometa", "**", "*.py"), recursive=True):
    exclude_dirs = ["validation"]
    if "__" in fpath:  # do not include __init__ and __main__.py files
        continue
    if os.path.basename(os.path.dirname(fpath)) in exclude_dirs:
        continue
    tmp_fpath = get_argparse_block(fpath)
    rst_dirpath = path_dir_docs_scripts(fpath)
    write_run_argparse_output(tmp_fpath, rst_dirpath, fpath)
    design_main_file(fpath, rst_dirpath)

link_write_usage()
remove_temp_files(tmp_fpath)


def main(args):
    pass


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="""
        Parses the argparse block of all autometa scripts, copy them to a 
        new file and run --help from there
        """)
    args = parser.parse_args()
    main(args)
