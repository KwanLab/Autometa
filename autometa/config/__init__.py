#!/usr/bin/env python
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

Utility functions for handling Autometa and user configuration.
"""


import logging
import os

from argparse import Namespace

from configparser import ConfigParser
from configparser import ExtendedInterpolation

logger = logging.getLogger(__name__)


DEFAULT_FPATH = os.path.join(os.path.dirname(__file__), "default.config")
AUTOMETA_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


def get_config(fpath):
    """Load the config provided at `fpath`.

    Parameters
    ----------
    fpath : str
        </path/to/file.config>

    Returns
    -------
    config.ConfigParser
        interpolated config object parsed from `fpath`.

    Raises
    -------
    FileNotFoundError
        Provided `fpath` does not exist.

    """
    # COMBAK: Checkpoint config
    if not os.path.exists(fpath) or os.stat(fpath).st_size == 0:
        raise FileNotFoundError(fpath)
    # https://stackoverflow.com/a/53274707/13118765
    converters = {"list": lambda x: [val.strip() for val in x.split(",")]}
    config = ConfigParser(interpolation=ExtendedInterpolation(), converters=converters)
    with open(fpath) as fh:
        config.read_file(fh)
    return config


DEFAULT_CONFIG = get_config(fpath=DEFAULT_FPATH)


def put_config(config, out):
    """Writes `config` to `out` and updates checkpoints checksum.

    Parameters
    ----------
    config : config.ConfigParser
        configuration containing user provided parameters and files information.
    out : str
        </path/to/output/file.config>

    Returns
    -------
    NoneType

    """
    with open(out, "w") as fh:
        config.write(fh)
    # COMBAK: Checkpoint config


def update_config(fpath, section, option, value):
    """Update `fpath` in `section` for `option` with `value`.

    Parameters
    ----------
    fpath : str
        </path/to/file.config>
    section : str
        `section` header to update within `fpath`.
    option : str
        `option` to update within `section`.
    value : str
        `value` to update `option`.

    Returns
    -------
    NoneType

    """
    cfg = get_config(fpath)
    cfg.set(section, option, value)
    put_config(cfg, fpath)
    logger.debug(f"updated {fpath} [{section}] option: {option} : {value}")


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
    type_converter = {
        "workspace": str,
        "project": int,
        "kingdoms": list,
        "metagenome_num": int,
        "length_cutoff": float,
        "cov_from_spades": bool,
        "kmer_size": int,
        "kmer_multiprocess": bool,
        "kmer_normalize": bool,
        "do_pca": bool,
        "pca_dims": int,
        "embedding_method": str,
        "do_taxonomy": bool,
        "taxon_method": str,
        "reversed": bool,
        "completeness": float,
        "purity": float,
        "binning_method": str,
        "verbose": bool,
        "force": bool,
        "usepickle": bool,
        "parallel": bool,
        "cpus": int,
        "config": str,
        "resume": bool,
        "fwd_reads": list,
        "rev_reads": list,
        "se_reads": list,
    }
    if fpath and (not os.path.exists(fpath) or os.stat(fpath).st_size == 0):
        raise FileNotFoundError(fpath)
    config = get_config(fpath) if fpath else DEFAULT_CONFIG
    namespace = Namespace()
    for section in config.sections():
        if section not in namespace:
            namespace.__dict__[section] = Namespace()
        for key, value in config.items(section):
            key = key.replace("-", "_")
            if section not in {"parameters", "files"} or key == "metagenomes":
                namespace.__dict__[section].__dict__[key] = value
                continue
            if type_converter.get(key) is not None:
                if type_converter.get(key) is bool:
                    value = config.getboolean(section, key)
                elif type_converter.get(key) is int:
                    value = config.getint(section, key)
                elif type_converter.get(key) is float:
                    value = config.getfloat(section, key)
                elif type_converter.get(key) is list:
                    value = config.getlist(section, key)
            namespace.__dict__[section].__dict__[key] = value
    return namespace


def set_home_dir():
    """Set the `home_dir` in autometa's default configuration (default.config)
    based on autometa's current location. If the `home_dir` variable is already
    set, then this will be used as the `home_dir` location.

    Returns
    -------
    str
        </path/to/package/autometa>

    """
    cfg = get_config(DEFAULT_FPATH)
    home_dir = cfg.get("common", "home_dir")
    if not os.path.isdir(home_dir):
        logger.info(
            f"Updated {os.path.basename(DEFAULT_FPATH)} ([common],home_dir): {AUTOMETA_DIR}"
        )
        update_config(DEFAULT_FPATH, "common", "home_dir", AUTOMETA_DIR)
        home_dir = AUTOMETA_DIR
    return home_dir
