#!/usr/bin/env python3
"""
Checks dependencies for Autometa Version 2
"""


import os
import sys

from config import databases

_EXECUTABLES = [
    'diamond',
    'hmmsearch',
    'hmmpress',
    'hmmscan',
    'prodigal',
    # 'makeblastdb',
    # 'blastp',
    'bowtie2',
    # 'md5sum',
    # 'gzip',
    # 'tar',
    # 'gunzip',
    # 'wget',
]


def which(program):
    """ Finds the full path for an executable and checks that read permissions exist
    See: https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    Arguments:
        program: the program to check
    Returns:
        The path if it was valid or None if not
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def check_executables(verbose=False):
    """
    Checks executable dependencies necessary to run autometa
    """
    # print('Checking Executables:')
    not_found = []
    for executable in _EXECUTABLES:
        found = which(executable)
        if not found:
            not_found.append(executable)
        if verbose:
            print(f'{executable}\t\t: {found}')
    if not_found:
        print(f'ExecutablesNotFound: {executable}')
        for executable in not_found:
            print(executable)
        sys.exit(1)
    print('Executable Dependencies Satisfied')
    return None

def load_databases(config_file, verbose=False):
    if not databases.check(config_file, verbose):
        sys.exit(1)
    print('Database Dependencies Satisfied')
    return databases.load(config_file)
