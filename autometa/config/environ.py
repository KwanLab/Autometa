#!/usr/bin/env python
"""
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

Configuration handling for Autometa environment.
"""


import logging
import os
import sys
import subprocess

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
    """Short summary.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    return {exe:which(exe) for exe in EXECUTABLES}

def diamond():
    """Get diamond version.

    Returns
    -------
    str
        version of diamond
    """
    exe = which('diamond')
    proc = subprocess.Popen(
        [exe,'version'],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL)
    stdout, stderr = proc.communicate()
    # stdout = b'diamond version 0.9.24\n'
    return stdout.decode().split()[-1]

def hmmsearch():
    """Get hmmsearch version.

    Returns
    -------
    str
        version of hmmsearch
    """
    exe = which('hmmsearch')
    proc = subprocess.Popen(
        [exe,'-h'],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode().split('#')[2]
    # stdout = ' HMMER 3.2.1 (June 2018); http://hmmer.org/\n'
    return stdout.strip().split()[1]

def hmmpress():
    """Get hmmpress version.

    Returns
    -------
    str
        version of hmmpress
    """
    exe = which('hmmpress')
    proc = subprocess.Popen(
        [exe,'-h'],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode().split('#')[2]
    # stdout = ' HMMER 3.2.1 (June 2018); http://hmmer.org/\n'
    return stdout.strip().split()[1]

def hmmscan():
    """Get hmmscan version.

    Returns
    -------
    str
        version of hmmscan
    """
    exe = which('hmmscan')
    proc = subprocess.Popen(
        [exe,'-h'],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL)
    stdout, stderr = proc.communicate()
    stdout = stdout.decode().split('#')[2]
    # stdout = ' HMMER 3.2.1 (June 2018); http://hmmer.org/\n'
    return stdout.strip().split()[1]

def prodigal():
    """Get prodigal version.

    Returns
    -------
    str
        version of prodigal
    """
    exe = which('prodigal')
    proc = subprocess.Popen(
        [exe,'-v'],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    # stderr = b'\nProdigal V2.6.3: February, 2016\n\n'
    return stderr.decode().strip().split(':')[0].replace('Prodigal V','')

def bowtie2():
    """Get bowtie2 version.

    Returns
    -------
    str
        version of bowtie2
    """
    exe = which('bowtie2')
    proc = subprocess.Popen(
        [exe,'--version'],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL)
    stdout, stderr = proc.communicate()
    # stdout 'bowtie2-align-s version 2.3.5\n64-bit\n
    return stdout.decode().split()[2]

def samtools():
    """Get samtools version.

    Returns
    -------
    str
        version of samtools
    """
    exe = which('samtools')
    proc = subprocess.Popen([exe], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    stderr = stderr.decode().strip().split('\n')[1]
    # stderr = 'Version: 1.10 (using htslib 1.10.2)'
    return stderr.split()[1]

def bedtools():
    """Get bedtools version.

    Returns
    -------
    str
        version of bedtools
    """
    exe = which('bedtools')
    proc = subprocess.Popen(
        [exe, '--version'],
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL)
    stdout, stderr = proc.communicate()
    # stdout = b'bedtools v2.29.2\n'
    return stdout.decode().strip().split()[-1].strip('v')

def get_versions(program=None):
    """Retrieve versions from all required executable dependencies.
    If `program` is provided will only return version for `program`.

    Parameters
    ----------
    program : str
        the program to retrieve the version.

    Returns
    -------
    dict
        {program:version, ...}

    Raises
    -------
    ValueError
        `program` is not a string
    KeyError
        `program` is not an executable dependency.

    """
    dispatcher = {
        'prodigal':prodigal,
        'diamond':diamond,
        'hmmsearch':hmmsearch,
        'hmmpress':hmmpress,
        'hmmscan':hmmscan,
        'prodigal':prodigal,
        'bowtie2':bowtie2,
        'samtools':samtools,
        'bedtools':bedtools,
    }
    if program:
        if type(program) is not str:
            raise ValueError(f'program is not string. given:{type(program)}')
        if program not in dispatcher:
            raise KeyError(f'{program} not in executables')
        return {program:dispatcher[program]()}
    versions = {}
    executables = find_executables()
    for exe,found in executables.items():
        if found:
            version = dispatcher[exe]()
        else:
            logger.warning(f'VersionUnavailable {exe}')
            version = 'ExecutableNotFound'
        versions.update({exe:version})
    return versions

def configure(config=DEFAULT_CONFIG):
    """Checks executable dependencies necessary to run autometa.
    Will update `config` with executable dependencies with details:
    1. presence/absence of dependency and its location
    2. versions

    Parameters
    ----------
    config : configparser.ConfigParser
        Description of parameter `config`.

    Returns
    -------
    2-tuple
        (config, satisfied)
        config updated with executables details
        Details:
        1. location of executable
        2. version of executable
        config : configparser.ConfigParser
        satisfied : bool

    """
    if not config.has_section('environ'):
        config.add_section('environ')
    if not config.has_section('versions'):
        config.add_section('versions')
    executables = find_executables()
    versions = get_versions()
    satisfied = True
    for executable,found in executables.items():
        version = versions.get(executable)
        if not config.has_option('environ', executable) and not found:
            satisfied = False
            logger.warning(f'executable not found: {executable}')
        elif not config.has_option('environ', executable):
            logger.debug(f'Updated executable {executable}: {found}')
            config.set('environ', executable, found)
            config.set('versions', executable, version)
        user_executable = config.get('environ', executable)
        if not which(user_executable):
            logger.debug(f'Updated executable {executable}: {found}')
            config.set('environ', executable, found)
            config.set('versions', executable, version)
        else:
            version = get_versions(user_executable).get(user_executable)
            config.set('versions', user_executable, version)
    return config, satisfied

def main(args):
    config,satisfied = configure(infpath=args.infpath)
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
