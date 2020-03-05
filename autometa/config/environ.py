#!/usr/bin/env python
"""
Configuration handling for Autometa environment.
"""


import logging
import os
import sys

from configparser import ConfigParser
from configparser import ExtendedInterpolation

from autometa.config import DEFAULT_CONFIG
from autometa.config import get_config
from autometa.config import put_config


logger = logging.getLogger(__name__)

EXECUTABLES = [
    'diamond',
    'hmmsearch',
    'hmmpress',
    'hmmscan',
    'prodigal',
    'bowtie2',
    'samtools',
    'bedtools',
]


def which(program):
    """Finds the full path for an executable and checks read permissions exist.

    See: https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python

    Returns:
        The path if it was valid or None if not

    Parameters
    ----------
    program : str
        the program to check

    Returns
    -------
    str
        </path/to/executable/> or ''

    Raises
    -------
    ExceptionName
        Why the exception is raised.

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
    return ''

def find_executables():
    return {exe:which(exe) for exe in EXECUTABLES}

def update_config(config):
    executables = find_executables()
    if not config.has_section('environ'):
        config.add_section('environ')
    satisfied = True
    for executable,found in executables.items():
        if not config.has_option('environ', executable) and not found:
            satisfied = False
            logger.warning(f'executable not found: {executable}')
        elif not config.has_option('environ', executable):
            logger.debug(f'Updated executable: {executable} : {found}')
            config.set('environ', executable, found)
        user_executable = config.get('environ', executable)
        if not which(user_executable):
            logger.debug(f'Updated executable: {executable} : {found}')
            config.set('environ', executable, found)
    logger.debug(f'Executable dependencies satisfied : {satisfied}')
    return config

def configure(config=DEFAULT_CONFIG):
    """
    Checks executable dependencies necessary to run autometa
    """
    return update_config(config)

def main(args):
    config = configure(infpath=args.infpath)
    if not args.out:
        import sys;sys.exit(0)
    put_config(config, args.out)

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Configure executables.config')
    parser.add_argument('--infpath',
        help='</path/to/output/executables.config>', required=True)
    parser.add_argument('--out',
        help='</path/to/output/executables.config>')
    args = parser.parse_args()
    main(args)
