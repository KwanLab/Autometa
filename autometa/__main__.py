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

from .config.user import AutometaUser


logger = logging.getLogger('autometa')


def init_logger(fpath=None, level=logging.INFO):
    """Initialize logger.

    By default will initialize streaming logger with DEBUG level messages.
    If `fpath` is provided, will write DEBUG level messages to `fpath` and
    set streaming messages to INFO.

    Parameters
    ----------
    fpath : str
        </path/to/file.log>
    level : int
        Overwrite default logging level behavior with provided `level`.
        This must be a constant from logging levels.
        See https://docs.python.org/3/library/logging.html#levels for details.
        i.e. logging.DEBUG, logging.INFO, etc. translates to 0,10, etc...

    Returns
    -------
    logging.Logger
        logging's Logger object to emit messages via methods:
        'warn','info','debug','error','exception','critical','fatal'

    Raises
    -------
    TypeError
        `level` must be an int
    ValueError
        `level` must be one of 0, 10, 20, 30, 40, 50
    """
    levels = {
        logging.NOTSET,
        logging.DEBUG,
        logging.INFO,
        logging.WARNING,
        logging.ERROR,
        logging.CRITICAL}
    if type(level) is not int:
        raise TypeError(f'{level} must be an int! {type(level)}')
    if level and level not in levels:
        raise ValueError(f'{level} not in levels: {levels}!')
    formatter = logging.Formatter(
        fmt='[%(asctime)s %(levelname)s] %(name)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')
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
    # Setup logger
    # timestamp = time.strftime("%Y-%m-%d_%H-%M-%S",time.gmtime())
    # log_fpath = args.log if args.log else f'{timestamp}_autometa.log'
    if args.debug:
        logger = init_logger(fpath=args.log, level=logging.DEBUG)
    else:
        logger = init_logger(fpath=args.log)
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
    import argparse
    import time
    cpus = mp.cpu_count()
    parser = argparse.ArgumentParser(description='Main script to run the Autometa pipeline.')
    parser.add_argument('config',
        help='Path to your metagenome.config file',
        nargs='*')
    parser.add_argument('--cpus',
        help=f'Num. cpus to use when updating/constructing databases (default: {cpus} cpus)',
        type=int,
        default=cpus)
    parser.add_argument('--debug',
        help='Stream debugging information to terminal',
        action='store_const',
        const=logging.DEBUG)
    parser.add_argument('--log', help='Path to write a log file (e.g. </path/to/autometa.log>)', type=str)
    parser.add_argument('--check-dependencies',
        help='Check user executables and databases accessible to Autometa and exit.',
        action='store_true')
    args = parser.parse_args()

    try:
        main(args)
    except KeyboardInterrupt:
        logger.info('User cancelled run. Exiting...')
    except Exception as err:
        issue_request = '''

        Please help us fix your problem!

        You may file an issue with us at https://github.com/KwanLab/Autometa/issues/new
        '''
        err.issue_request = issue_request
        logger.exception(err)
        logger.info(err.issue_request)

if __name__ == '__main__':
    entrypoint()
