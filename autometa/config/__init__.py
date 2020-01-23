#!/usr/bin/env python

import os

from argparse import Namespace


from configparser import ConfigParser
from configparser import ExtendedInterpolation

DEFAULT_FPATH = os.path.join(os.path.dirname(__file__), 'default.config')
AUTOMETA_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
PROJECTS_DIR = os.path.join(AUTOMETA_DIR, 'projects')



def get_config(fpath):
    if not os.path.exists(fpath) or os.stat(fpath).st_size == 0:
        raise FileNotFoundError(fpath)
    config = ConfigParser(interpolation=ExtendedInterpolation())
    with open(fpath) as fh:
        config.read_file(fh)
    return config

def put_config(config, out):
    with open(out, 'w') as fh:
        config.write(fh)

def update_config(fpath, section, option, value):
    c = get_config(fpath)
    c.set(section,option,value)
    put_config(c, fpath)
    logger.debug(f'updated {fpath} [{section}] option: {option} : {value}')

DEFAULT_CONFIG = get_config(fpath=DEFAULT_FPATH)

parameters = {
     'projects':str,
     'project':int,
     'add-metagenome':str,
     'resume':int,
     'length_cutoff':float,
     'cov_from_spades':bool,
     'kmer_size':int,
     'kmer_multiprocess':bool,
     'kmer_normalize':bool,
     'do_pca':bool,
     'pca_dims':int,
     'embedding_method':str,
     'taxon_method':str,
     'reversed':bool,
     'completeness':float,
     'purity':float,
     'binning_method':str,
     'verbose':bool,
     'force':bool,
     'usepickle':bool,
     'parallel':bool,
     'cpus':int,
     'config':str,
     'resume':bool,
}


def parse_config(fpath=None):
    """Generate argparse namespace (args) from config file.

    Parameters
    ----------
    fpath : str
        </path/to/file.config> (default is DEFAULT_CONFIG in autometa.config)

    Returns
    -------
    argparse.Namespace
        namespace typical to parser.parse_args() method from argparse

    Raises
    -------
    FileNotFoundError
        provided `fpath` does not exist.

    """
    if fpath and (not os.path.exists(fpath) or os.stat(fpath).st_size == 0):
        raise FileNotFoundError(fpath)
    config = get_config(fpath) if fpath else DEFAULT_CONFIG
    namespace = Namespace()
    for section in config.sections():
        if section not in namespace:
            namespace.__dict__[section] = Namespace()
        for key, value in config.items(section):
            key = key.replace('-', '_')
            if section != 'parameters' or key == 'metagenomes':
                namespace.__dict__[section].__dict__[key] = value
                continue
            if parameters.get(key) is not None:
                if parameters.get(key) is bool:
                    value = config.getboolean(section,key)
                elif parameters.get(key) is int:
                    value = config.getint(section, key)
                elif parameters.get(key) is float:
                    value = config.getfloat(section,key)
            namespace.__dict__[section].__dict__[key] = value
    return namespace
