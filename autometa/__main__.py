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

Main script to run Autometa
"""


import logging
import os
import sys

import multiprocessing as mp

from autometa.config.user import AutometaUser


logger = logging.getLogger("autometa")


def init_logger(fpath=None, verbosity=0):
    """Initialize logger.

    By default will initialize streaming logger with INFO level messages.
    If `fpath` is provided, will write DEBUG level messages to `fpath` and
    set streaming messages to INFO.

    Parameters
    ----------
    fpath : str, optional
        </path/to/file.log>
    verbosity : int, optional
        Overwrite default logging level behavior with provided `verbosity`.
        This must be between 0-2 where 0 is the least information and 2 is the most information.
        See https://docs.python.org/3/library/logging.html#levels for details.

    Returns
    -------
    logging.Logger
        logging's Logger object to emit messages via methods:
        'warn','info','debug','error','exception','critical','fatal'

    Raises
    -------
    TypeError
        `verbosity` must be an int
    ValueError
        `verbosity` must be between of 0 and 2
    """
    log_levels = {
        0: logging.WARN,
        1: logging.INFO,
        2: logging.DEBUG,
    }

    if type(verbosity) is not int:
        raise TypeError(f"{verbosity} must be an int! {type(verbosity)}")
    if verbosity and verbosity not in log_levels:
        raise ValueError(f"{verbosity} not in log_levels: {log_levels}!")

    level = log_levels.get(verbosity)

    formatter = logging.Formatter(
        fmt="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    # Construct file/stream logging handlers
    streamhandler = logging.StreamHandler()
    streamhandler.setFormatter(formatter)
    if fpath:
        filehandler = logging.FileHandler(fpath)
        filehandler.setFormatter(formatter)
        logger.addHandler(filehandler)

    streamhandler.setLevel(level)
    logger.addHandler(streamhandler)
    logger.setLevel(logging.DEBUG)
    return logger


def main(args):
    """
    Main logic for running autometa pipeline.

    Warning: This should be called by `entrypoint` and not directly.

    Parameters
    ----------
    args : argparse.Namespace
        namespace containing config information for AutometaUser

    Returns
    -------
    NoneType
        Nothing if no errors are encountered.

    """

    logger = init_logger(fpath=args.log, verbosity=args.verbosity)
    # Configure AutometaUser
    # TODO: master from WorkQueue is AutometaUser
    user = AutometaUser(nproc=args.cpus)
    user.configure(dryrun=args.check_dependencies)

    for config in args.config:
        # TODO: Add directions to master from WorkQueue
        mgargs = user.prepare_binning_args(config)
        user.run_binning(mgargs)
        # user.refine_binning()
        # user.process_binning()
    # user.get_pangenomes()


def entrypoint():
    """
    Main entrypoint for autometa pipeline.

    Note, a requirement of packaging and distribution is for entrypoints to not
    require any arguments. This is a wrapper to the main functionality of running
    autometa via the `main` function.

    Returns
    -------
    NoneType

    """
    import argparse

    # import time
    cpus = mp.cpu_count()
    parser = argparse.ArgumentParser(
        description="Main script to run the Autometa pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("config", help="Path to your metagenome.config file", nargs="*")
    parser.add_argument(
        "--cpus",
        help="Num. cpus to use when updating/constructing databases",
        type=int,
        default=cpus,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbosity",
        action="count",
        default=0,
        help="Verbosity (between 1-2 occurrences with more leading to more "
        "verbose logging). WARN=0, INFO=1, DEBUG=2",
    )
    parser.add_argument(
        "--log",
        help="Path to write a log file (e.g. </path/to/autometa.log>)",
        type=str,
    )
    parser.add_argument(
        "--check-dependencies",
        help="Check user executables and databases accessible to Autometa and exit.",
        action="store_true",
    )
    args = parser.parse_args()

    try:
        main(args)
    except KeyboardInterrupt:
        logger.info("User cancelled run. Exiting...")
        sys.exit(1)
    except Exception as err:
        issue_request = """
        An error was encountered!

        Please help us fix your problem!

        You may file an issue with us at https://github.com/KwanLab/Autometa/issues/new/choose
        """
        err.issue_request = issue_request
        logger.exception(err)
        logger.info(issue_request)
        sys.exit(1)


if __name__ == "__main__":
    entrypoint()
