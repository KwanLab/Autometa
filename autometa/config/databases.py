#!/usr/bin/env python3
"""
Database configuration handling for autometa
"""
from configparser import ConfigParser
from configparser import ExtendedInterpolation
import os

_DEFAULT_NAME = 'default.cfg'
_BASEDIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG = os.path.join(_BASEDIR, _DEFAULT_NAME)
_NCBI_OPTS = [
    'nodes',
    'names',
    'accession2taxid',
    'blastdb',
]
_MARKER_OPTS = [
    'bacteria_single_copy',
    'bacteria_single_copy_cutoffs',
    'archaea_single_copy',
    'archaea_single_copy_cutoffs',
]
_DATABASES = {
    'ncbi':_NCBI_OPTS,
    'markers':_MARKER_OPTS,
}


def check(config_fpath=DEFAULT_CONFIG, verbose=False):
    passed = True
    # load generic configuration settings
    config = ConfigParser(interpolation=ExtendedInterpolation())
    with open(config_fpath, 'r') as handle:
        config.read_file(handle)
    for section,options in _DATABASES.items():
        has_section = section in config.sections()
        if verbose:
            print('{} section exists: {}'.format(section, has_section))
        if not has_section:
            passed = False
        for candidate in options:
            has_option = config.has_option(section, candidate)
            if has_option:
                fpath = config.get(section,candidate)
                exists = os.path.exists(fpath)
            else:
                False
            if verbose:
                print(
                    f'{section} has option {candidate}:{has_option}'
                    f' FileExists:{exists}'
                )
            if not has_option or not exists:
                passed = False
    return passed

def load(config_fpath=DEFAULT_CONFIG):
    # load generic configuration settings
    config = ConfigParser(interpolation=ExtendedInterpolation())
    with open(config_fpath, 'r') as handle:
        config.read_file(handle)
    return config
