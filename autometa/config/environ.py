#!/usr/bin/env python
"""
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

Configuration handling for Autometa environment.
"""


import logging
import os
import shutil
import subprocess

from configparser import ConfigParser
from typing import Dict, Tuple, Union


logger = logging.getLogger(__name__)

EXECUTABLES = [
    "diamond",
    "hmmsearch",
    "hmmpress",
    "hmmscan",
    "prodigal",
    "bowtie2",
    "samtools",
    "bedtools",
]


def find_executables():
    """Retrieves executable file paths by looking in Autometa dependent executables.

    Returns
    -------
    dict
        {executable:</path/to/executable>, ...}

    """
    return {exe: shutil.which(exe, mode=os.X_OK) for exe in EXECUTABLES}


def diamond():
    """Get diamond version.

    Returns
    -------
    str
        version of diamond

    """
    exe = shutil.which("diamond", mode=os.X_OK)
    proc = subprocess.Popen(
        [exe, "version"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
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
    exe = shutil.which("hmmsearch", mode=os.X_OK)
    proc = subprocess.Popen(
        [exe, "-h"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    stdout, stderr = proc.communicate()
    stdout = stdout.decode().split("#")[2]
    # stdout = ' HMMER 3.2.1 (June 2018); http://hmmer.org/\n'
    return stdout.strip().split()[1]


def hmmpress():
    """Get hmmpress version.

    Returns
    -------
    str
        version of hmmpress
    """
    exe = shutil.which("hmmpress", mode=os.X_OK)
    proc = subprocess.Popen(
        [exe, "-h"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    stdout, stderr = proc.communicate()
    stdout = stdout.decode().split("#")[2]
    # stdout = ' HMMER 3.2.1 (June 2018); http://hmmer.org/\n'
    return stdout.strip().split()[1]


def hmmscan():
    """Get hmmscan version.

    Returns
    -------
    str
        version of hmmscan
    """
    exe = shutil.which("hmmscan", mode=os.X_OK)
    proc = subprocess.Popen(
        [exe, "-h"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    stdout, stderr = proc.communicate()
    stdout = stdout.decode().split("#")[2]
    # stdout = ' HMMER 3.2.1 (June 2018); http://hmmer.org/\n'
    return stdout.strip().split()[1]


def prodigal():
    """Get prodigal version.

    Returns
    -------
    str
        version of prodigal
    """
    exe = shutil.which("prodigal", mode=os.X_OK)
    proc = subprocess.Popen(
        [exe, "-v"], stdout=subprocess.DEVNULL, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    # stderr = b'\nProdigal V2.6.3: February, 2016\n\n'
    return stderr.decode().strip().split(":")[0].replace("Prodigal V", "")


def bowtie2():
    """Get bowtie2 version.

    Returns
    -------
    str
        version of bowtie2
    """
    exe = shutil.which("bowtie2", mode=os.X_OK)
    proc = subprocess.Popen(
        [exe, "--version"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
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
    exe = shutil.which("samtools", mode=os.X_OK)
    proc = subprocess.Popen([exe], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    stderr = stderr.decode().strip().split("\n")[1]
    # stderr = 'Version: 1.10 (using htslib 1.10.2)'
    return stderr.split()[1]


def bedtools():
    """Get bedtools version.

    Returns
    -------
    str
        version of bedtools
    """
    exe = shutil.which("bedtools", mode=os.X_OK)
    proc = subprocess.Popen(
        [exe, "--version"], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL
    )
    stdout, stderr = proc.communicate()
    # stdout = b'bedtools v2.29.2\n'
    return stdout.decode().strip().split()[-1].strip("v")


def get_versions(program: str = None) -> Union[Dict[str, str], str]:
    """
    Retrieve versions from all required executable dependencies.
    If `program` is provided will only return version for `program`.

    See: https://stackoverflow.com/a/834451/12671809

    Parameters
    ----------
    program : str, optional
        the program to retrieve the version, by default None

    Returns
    -------
    dict or str
        if program is None: dict - {program:version, ...}
        if program: str - version

    Raises
    -------
    ValueError
        `program` is not a string
    KeyError
        `program` is not an executable dependency.
    """
    if program:
        if not isinstance(program, str):
            raise TypeError(f"program is not string. given:{type(program)}")
        exe_name = os.path.basename(program)
        if exe_name not in globals():
            raise KeyError(f"{exe_name} not in executables")
        try:
            return globals()[exe_name]()
        except TypeError:
            logger.warning(
                f"{exe_name} not found. This may impact a stage of the Autometa pipeline."
            )
            return "Not found"
    versions = {}
    executables = find_executables()
    for exe, found in executables.items():
        if found:
            # get_version accesses the location of the function for that program, eg. find the location of bedtools()
            # get_version() wraps the location and execute that function, eg. execute bedtools()
            get_version = globals()[exe]
            version = get_version()
        else:
            logger.warning(f"VersionUnavailable {exe}")
            version = "ExecutableNotFound"
        versions.update({exe: version})
    return versions


def configure(config: ConfigParser) -> Tuple[ConfigParser, bool]:
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
    if not config.has_section("environ"):
        config.add_section("environ")
    if not config.has_section("versions"):
        config.add_section("versions")
    executables = find_executables()
    versions = get_versions()
    satisfied = True
    for executable, found in executables.items():
        version = versions.get(executable)
        if not config.has_option("environ", executable) and not found:
            satisfied = False
            logger.warning(f"executable not found: {executable}")
        elif not config.has_option("environ", executable):
            logger.debug(f"{executable}: {found} (version: {version})")
            config.set("environ", executable, found)
            config.set("versions", executable, version)
        user_executable = config.get("environ", executable)
        if not shutil.which(user_executable, mode=os.X_OK):
            logger.debug(f"{executable}: {found} (version: {version})")
            config.set("environ", executable, str(found))
            config.set("versions", executable, version)
        else:
            version = get_versions(user_executable)
            if not os.path.basename(user_executable) in config.options("versions"):
                config.set("versions", user_executable, version)
    return config, satisfied


if __name__ == "__main__":
    import sys

    sys.exit(0)
