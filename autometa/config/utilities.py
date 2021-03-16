#!/usr/bin/env python


import os
import logging

from argparse import Namespace
from configparser import ConfigParser
from configparser import ExtendedInterpolation
from autometa.config import environ

logger = logging.getLogger(__name__)


DEFAULT_FPATH = os.path.join(os.path.dirname(__file__), "default.config")
AUTOMETA_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))


def get_config(fpath: str) -> ConfigParser:
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
    if not os.path.exists(fpath) or not os.path.getsize(fpath):
        raise FileNotFoundError(fpath)
    # https://stackoverflow.com/a/53274707/13118765
    converters = {"list": lambda x: [val.strip() for val in x.split(",")]}
    config = ConfigParser(interpolation=ExtendedInterpolation(), converters=converters)
    with open(fpath) as fh:
        config.read_file(fh)
    return config


DEFAULT_CONFIG = get_config(fpath=DEFAULT_FPATH)


def put_config(config: ConfigParser, out: str) -> None:
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


def update_config(
    section: str, option: str, value: str, fpath: str = DEFAULT_FPATH
) -> None:
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
    if not cfg.has_option(section, option):
        options = ", ".join(cfg.options(section))
        raise KeyError(
            f"option: `{option}` not available in section: `{section}`. Available options: {options}"
        )
    cfg.set(section, option, value)
    put_config(cfg, fpath)
    logger.debug(f"updated {fpath} [{section}] option: {option} : {value}")


def parse_args(fpath: str = None) -> Namespace:
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
        "kmer_transform": str,
        "do_pca": bool,
        "pca_dims": int,
        "embed_dimensions": int,
        "embed_method": str,
        "do_taxonomy": bool,
        "tmpdir": str,
        "taxon_method": str,
        "reverse_ranks": bool,
        "starting_rank": str,
        "clustering_method": str,
        "completeness": float,
        "purity": float,
        "binning_method": str,
        "verbose": bool,
        "force": bool,
        "usepickle": bool,
        "parallel": bool,
        "cpus": int,
        "seed": int,
        "config": str,
        "resume": bool,
        "fwd_reads": list,
        "rev_reads": list,
        "se_reads": list,
    }
    if fpath and (not os.path.exists(fpath) or not os.path.getsize(fpath)):
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


def set_home_dir() -> str:
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
            f"Updating {os.path.basename(DEFAULT_FPATH)} ([common],home_dir): {AUTOMETA_DIR}"
        )
        update_config(section="common", option="home_dir", value=AUTOMETA_DIR)
        home_dir = AUTOMETA_DIR
    return home_dir


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="Update Autometa configuration using provided arguments"
    )
    logging_group = parser.add_argument_group("Logging")
    update_group = parser.add_argument_group("Updating")
    update_group.add_argument(
        "--section",
        help="config section to update",
        choices=["environ", "databases", "ncbi", "markers"],
        required=False,
    )
    update_group.add_argument(
        "--option", help="option in `--section` to update", required=False
    )
    update_group.add_argument(
        "--value", help="Value to update `--option`", required=False
    )
    logging_group.add_argument(
        "--print", help="Print configuration without updating", action="store_true"
    )
    args = parser.parse_args()

    # First update home_dir option in common section with Autometa installation path
    set_home_dir()
    # Now ensure executables are available: save versions to config
    cfg = get_config(DEFAULT_FPATH)
    cfg, environ_satisfied = environ.configure(cfg)
    put_config(cfg, DEFAULT_FPATH)
    logger.debug(f"environment dependencies satisifed: {environ_satisfied}")
    if args.print:
        print("section\toption\tvalue")
        for section in cfg.sections():
            for option in cfg.options(section):
                value = cfg.get(section, option)
                print(f"{section}\t{option}\t{value}")
        exit(0)
    for option in [args.section, args.option, args.value]:
        if not option:
            raise ValueError(
                "the following arguments are required: --section, --option, --value"
            )
    update_config(section=args.section, option=args.option, value=args.value)


if __name__ == "__main__":
    main()
