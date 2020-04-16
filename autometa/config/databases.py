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

Configuration handling for Autometa Databases.
"""


import logging
import os
import requests

import multiprocessing as mp

from configparser import ConfigParser
from configparser import ExtendedInterpolation
from ftplib import FTP

from autometa.config import get_config
from autometa.config import DEFAULT_FPATH
from autometa.config import DEFAULT_CONFIG
from autometa.config import put_config
from autometa.config import AUTOMETA_DIR
from autometa.common.utilities import untar
from autometa.common.external import diamond
from autometa.common.external import hmmer


logger = logging.getLogger(__name__)
urllib_logger = logging.getLogger("urllib3")
urllib_logger.setLevel(logging.WARNING)


class Databases:
    """docstring for Databases."""
    SECTIONS = {
        'ncbi':[
            'nodes',
            'names',
            'merged',
            'accession2taxid',
            'nr',
        ],
        'markers':[
            'bacteria_single_copy',
            'bacteria_single_copy_cutoffs',
            'archaea_single_copy',
            'archaea_single_copy_cutoffs',
        ],
    }
    def __init__(self, config=DEFAULT_CONFIG, dryrun=False, nproc=mp.cpu_count()):
        if type(config) is not ConfigParser:
            raise TypeError(f'config is not ConfigParser : {type(config)}')
        if type(dryrun) is not bool:
            raise TypeError(f'dryrun must be True or False. type: {type(dryrun)}')

        self.config=config
        self.dryrun=dryrun
        self.nproc=nproc
        self.prepare_sections()
        self.ncbi_dir = self.config.get('databases','ncbi')
        self.markers_dir = self.config.get('databases','markers')

    @property
    def satisfied(self):
        return self.get_missing(validate=True)

    def prepare_sections(self):
        """Add database sections to 'databases' if missing.

        Returns
        -------
        NoneType

        """
        if not self.config.has_section('databases'):
            self.config.add_section('databases')
        for section in Databases.SECTIONS:
            if self.config.has_option('databases', section):
                continue
            outdir = DEFAULT_CONFIG.get('databases', section)
            self.config.set('databases', section, outdir)

    def format_nr(self):
        """Format NCBI nr.gz database into diamond formatted database nr.dmnd.

        Returns
        -------
        NoneType
            config updated option:'nr' in section:'ncbi'.

        """
        db_infpath = self.config.get('ncbi','nr')
        db_outfpath = db_infpath.replace('.gz','.dmnd')
        if not self.dryrun and not os.path.exists(db_outfpath):
            diamond.makedatabase(fasta=nr, database=db_infpath, nproc=self.nproc)
        self.config.set('ncbi','nr', db_outfpath)
        logger.debug(f'set ncbi nr: {db_outfpath}')

    def extract_taxdump(self):
        """Extract autometa required files from ncbi taxdump directory.

        Extracts nodes.dmp, names.dmp and merged.dmp from taxdump.tar.gz

        Returns
        -------
        NoneType

        """
        taxdump = self.config.get('ncbi','taxdump')
        taxdump_files = [
            ('nodes','nodes.dmp'),
            ('names','names.dmp'),
            ('merged','merged.dmp'),
        ]
        for option,fname in taxdump_files:
            outfpath = os.path.join(self.ncbi_dir,fname)
            if not self.dryrun and not os.path.exists(outfpath):
                outfpath = untar(taxdump, self.ncbi_dir, fname)
            logger.debug(f'UPDATE (ncbi,{option}): {outfpath}')
            self.config.set('ncbi',option,outfpath)

    def update_ncbi(self, options):
        """Update NCBI database files (taxdump.tar.gz and nr.gz).

        Parameters
        ----------
        options : set
            Set of options to update

        Returns
        -------
        NoneType
            Description of returned object.

        Raises
        -------
        ConnectionError
            NCBI file download failed.

        """
        # Download required NCBI database files
        if not os.path.exists(self.ncbi_dir):
            os.makedirs(self.ncbi_dir)
        host = DEFAULT_CONFIG.get('ncbi','host')
        for option in options:
            ftp_fullpath = DEFAULT_CONFIG.get('database_urls',option)
            ftp_fpath = ftp_fullpath.split(host)[-1]
            if self.config.has_option('ncbi', option):
                outfpath = self.config.get('ncbi', option)
            else:
                outfname = os.path.basename(ftp_fpath)
                outfpath = os.path.join(self.ncbi_dir, outfname)
            if not self.dryrun and not os.path.exists(outfpath):
                with FTP(host) as ftp, open(outfpath, 'wb') as fp:
                    ftp.login()
                    logger.debug(f'starting {option} download')
                    result = ftp.retrbinary(f'RETR {ftp_fpath}', fp.write)
                    if not result.startswith('226 Transfer complete'):
                        raise ConnectionError(f'{option} download failed')
                    ftp.quit()
            logger.debug(f'UPDATE: (ncbi,{option}): {outfpath}')
            self.config.set('ncbi', option, outfpath)
        # Extract/format respective NCBI files
        self.extract_taxdump()
        self.format_nr()

    def update_markers(self, options):
        """Update single-copy markers hmms and cutoffs.

        Parameters
        ----------
        options : set
            Description of parameter `options`.

        Returns
        -------
        NoneType

        Raises
        -------
        ConnectionError
            marker file download failed.

        """
        if not os.path.exists(self.markers_dir):
            os.makedirs(self.markers_dir)
        for option in options:
            url = DEFAULT_CONFIG.get('database_urls', option)
            if self.config.has_option('markers', option):
                outfpath = self.config.get('markers', option)
            else:
                outfname = os.path.basename(url)
                outfpath = os.path.join(self.markers_dir, outfname)
            if self.dryrun:
                logger.debug(f'UPDATE: (markers,{option}): {outfpath}')
                self.config.set('markers', option, outfpath)
                continue
            with requests.Session() as session:
                resp = session.get(url)
            if not resp.ok:
                raise ConnectionError(f'Failed to retrieve {url}')
            with open(outfpath, 'w') as outfh:
                outfh.write(resp.text)
            self.config.set('markers', option, outfpath)
            if outfpath.endswith('.hmm'):
                hmmer.hmmpress(outfpath)

    def get_missing(self, validate=False):
        """Retrieve all database files from all database sections that are not
        available.

        Returns
        -------
        bool or dict

            - if `validate` is True : bool

                all available evaluates to True, otherwise False

            - if `validate` is False : dict

                {section:{option, option,...}, section:{...}, ...}

        """
        missing = {}
        for section in Databases.SECTIONS:
            for option in self.config.options(section):
                if option not in Databases.SECTIONS.get(section):
                    # Skip user added options not required by Autometa
                    continue
                fpath = self.config.get(section,option)
                if os.path.exists(fpath) and os.stat(fpath).st_size >= 0:
                    # TODO: [Checkpoint validation]
                    logger.debug(f'({section},{option}): {fpath}')
                    continue
                if validate:
                    return False
                if section in missing:
                    missing[section].add(option)
                else:
                    missing.update({section:set([option])})
        for section,opts in missing.items():
            for opt in opts:
                logger.debug(f'MISSING: ({section},{opt})')
        return True if validate else missing

    def update_missing(self):
        """Download and format databases for all options in each section.

        NOTE: This will only perform the download and formatting if self.dryrun is False

        Returns
        -------
        NoneType
            config updated with required missing sections.

        """
        dispatcher = {'ncbi':self.update_ncbi, 'markers':self.update_markers}
        missing = self.get_missing()
        for section,options in missing.items():
            if section == 'ncbi':
                if 'nodes' in options or 'names' in options or 'merged' in options:
                    options.discard('nodes')
                    options.discard('names')
                    options.discard('merged')
                    options.add('taxdump')
            dispatcher[section](options)

    def configure(self):
        """Checks database files

        Returns
        -------
        configparser.ConfigParser
            config with updated options in respective databases sections.

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        self.update_missing()
        return self.config

def main():
    import argparse
    import logging as logger

    cpus = mp.cpu_count()
    logger.basicConfig(
        format='[%(asctime)s %(levelname)s] %(name)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)

    parser = argparse.ArgumentParser('databases config', epilog='By default, with no arguments, will download/format databases into default databases directory.')
    parser.add_argument('--config', help='</path/to/input/database.config>', default=DEFAULT_FPATH)
    parser.add_argument('--dryrun', help='Log configuration actions but do not perform them.',
        action='store_true', default=False)
    parser.add_argument('--nproc',
        help=f'num. cpus to use for DB formatting. (default {cpus})', type=int, default=cpus)
    parser.add_argument('--out', help='</path/to/output/database.config>')
    args = parser.parse_args()

    config = get_config(args.config)
    dbs = Databases(config=config, dryrun=args.dryrun, nproc=args.nproc)
    logger.debug(f'Configuring databases')
    config = dbs.configure()
    dbs = Databases(config=config, dryrun=args.dryrun, nproc=args.nproc)
    logger.info(f'Database dependencies satisfied: {dbs.satisfied}')
    if not args.out:
        import sys;sys.exit(0)
    put_config(config, args.out)
    logger.debug(f'{args.out} written.')

if __name__ == '__main__':
    main()
