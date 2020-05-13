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

This file contains the Databases class responsible for configuration handling
of Autometa Databases.
"""


import logging
import os
import requests
import tempfile

import multiprocessing as mp

from configparser import ConfigParser
from configparser import ExtendedInterpolation
from ftplib import FTP
from glob import glob

from autometa.config import get_config
from autometa.config import DEFAULT_FPATH
from autometa.config import DEFAULT_CONFIG
from autometa.config import put_config
from autometa.config import AUTOMETA_DIR
from autometa.common.utilities import untar
from autometa.common.utilities import calc_checksum
from autometa.common.utilities import read_checksum
from autometa.common.utilities import write_checksum
from autometa.common.external import diamond
from autometa.common.external import hmmer


logger = logging.getLogger(__name__)
urllib_logger = logging.getLogger("urllib3")
urllib_logger.setLevel(logging.WARNING)


class Databases:
    """Database class containing methods to allow downloading/formatting/updating
    Autometa database dependencies.

    Parameters
    ----------
    config : config.ConfigParser
        Config containing database dependency information (the default is DEFAULT_CONFIG).
    dryrun : bool
        Run through database checking without performing downloads/formatting (the default is False).
    nproc : int
        Number of processors to use to perform database formatting. (the default is mp.cpu_count()).
    do_update : bool
        Overwrite existing databases with more up-to-date database files. (the default is False).

    Attributes
    ----------
    ncbi_dir : str
        </path/to/databases/ncbi>
    markers_dir : str
        </path/to/databases/markers>
    SECTIONS : dict
        keys are `sections` respective to database config sections and values
        are options within the `sections`.

    """

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
    def __init__(self, config=DEFAULT_CONFIG, dryrun=False, nproc=mp.cpu_count(),
        do_update=False):
        if not isinstance(config, ConfigParser):
            raise TypeError(f'config is not ConfigParser : {type(config)}')
        if not isinstance(dryrun, bool) and isins:
            raise TypeError(f'dryrun must be boolean. type: {type(dryrun)}')

        self.config=config
        self.dryrun=dryrun
        self.nproc=nproc
        self.do_update = do_update
        if not self.config.has_section('databases'):
            self.config.add_section('databases')
        for section in Databases.SECTIONS:
            if self.config.has_option('databases', section):
                continue
            outdir = DEFAULT_CONFIG.get('databases', section)
            self.config.set('databases', section, outdir)
        self.ncbi_dir = self.config.get('databases', 'ncbi')
        self.markers_dir = self.config.get('databases', 'markers')
        for outdir in {self.ncbi_dir, self.markers_dir}:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

    def satisfied(self, section=None, compare_checksums=False):
        """Determines whether all database dependencies are satisfied.

        Parameters
        ----------
        section : str
            section to retrieve for `checksums` section.
            Choices include: 'ncbi' and 'markers'.
        compare_checksums : bool, optional
            Also check if database information is up-to-date with current hosted databases.
            (default is False).

        Returns
        -------
        bool
            True if all database dependencies are satisfied, otherwise False.

        """
        any_missing = self.get_missing(section=section, validate=True)
        if compare_checksums:
            any_invalid = self.compare_checksums(section=section, validate=True)
        else:
            any_invalid = False
        if not any_missing and not any_invalid:
            return True
        else:
            return False

    def get_remote_checksum(self, section, option):
        """Get the checksum from provided `section` respective to `option` in `self.config`.

        Parameters
        ----------
        section : str
            section to retrieve for `checksums` section.
            Choices include: 'ncbi' and 'markers'.
        option : str
            `option` in `checksums` section corresponding to the section checksum file.

        Returns
        -------
        str
            checksum of remote md5 file. e.g. 'hash filename\n'

        Raises
        -------
        ConnectionError
            Failed to connect to host for provided `option`.

        """
        if section not in {'ncbi','markers'}:
            raise ValueError(f'"section" must be "ncbi" or "markers". Provided: {section}')
        if section == 'ncbi':
            host = self.config.get(section, 'host')
            ftp_fullpath = self.config.get('checksums', option)
            chksum_fpath = ftp_fullpath.split(host)[-1]
            with FTP(host) as ftp, tempfile.TemporaryFile() as fp:
                ftp.login()
                result = ftp.retrbinary(f'RETR {chksum_fpath}', fp.write)
                if not result.startswith('226 Transfer complete'):
                    raise ConnectionError(f'{chksum_fpath} download failed')
                ftp.quit()
                fp.seek(0)
                checksum = fp.read().decode()
        elif section == 'markers':
            url = self.config.get('checksums', option)
            with requests.Session() as session:
                resp = session.get(url)
                if not resp.ok:
                    raise ConnectionError(f'Failed to retrieve {url}')
                checksum = resp.text
        return checksum

    def format_nr(self):
        """Construct a diamond formatted database (nr.dmnd) from `nr` option
        in `ncbi` section in user config.

        NOTE: The checksum 'nr.gz.md5' will only be generated if nr.dmnd
        construction is successful. If the provided `nr` option in `ncbi` is 'nr.gz'
        the database will be removed after successful database formatting.

        Returns
        -------
        NoneType
            config updated option:'nr' in section:'ncbi'.

        """
        db_infpath = self.config.get('ncbi','nr')
        db_infpath_md5 = f'{db_infpath}.md5'
        db_outfpath = db_infpath.replace('.gz','.dmnd')
        checksums_match = False
        if os.path.exists(db_infpath_md5) and os.path.exists(db_outfpath):
            # If nr.dmnd.md5 exists, then we will check if nr.gz.md5 is matching
            current_checksum = read_checksum(db_infpath_md5)
            current_hash, __ = current_checksum.split()
            host = self.config.get('ncbi', 'host')
            remote_checksum = self.get_remote_checksum(host, 'nr')
            remote_hash, __ = remote_checksum.split()
            if current_hash == remote_hash:
                logger.debug(f'nr checksums match, skipping...')
                checksums_match = True

        no_update_db = os.path.exists(db_outfpath) and not self.do_update
        if self.dryrun or no_update_db or checksums_match:
            self.config.set('ncbi', 'nr', db_outfpath)
            logger.debug(f'set ncbi nr: {db_outfpath}')
            return

        diamond.makedatabase(fasta=db_infpath, database=db_outfpath, nproc=self.nproc)
        if os.path.basename(db_infpath) == 'nr.gz':
            # nr.gz will be removed after successful nr.dmnd construction
            write_checksum(db_infpath, db_infpath_md5)
            os.remove(db_infpath)

        self.config.set('ncbi', 'nr', db_outfpath)
        logger.debug(f'set ncbi nr: {db_outfpath}')

    def extract_taxdump(self):
        """Extract autometa required files from ncbi taxdump.tar.gz archive
        into ncbi databases directory and update user config with extracted paths.

        This only extracts nodes.dmp, names.dmp and merged.dmp from
        taxdump.tar.gz if the files do not already exist. If `do_update`
        was originally supplied as `True` to the Databases instance, then the
        previous files will be replaced by the new taxdump files.

        After successful extraction of the files, a checksum will be written of the
        archive for future checking and then the archive will be removed to save user disk space.

        Returns
        -------
        NoneType
            Will update `self.config` section `ncbi` with options 'nodes','names','merged'

        """
        taxdump_fpath = self.config.get('ncbi','taxdump')
        taxdump_files = [
            ('nodes','nodes.dmp'),
            ('names','names.dmp'),
            ('merged','merged.dmp'),
        ]
        for option,fname in taxdump_files:
            outfpath = os.path.join(self.ncbi_dir, fname)
            if self.dryrun:
                logger.debug(f'UPDATE (ncbi,{option}): {outfpath}')
                self.config.set('ncbi', option, outfpath)
                continue
            # Only update the taxdump files if the user says to do an update.
            if self.do_update and os.path.exists(outfpath):
                os.remove(outfpath)
            # Only extract the taxdump files if this is not a "dryrun"
            if not os.path.exists(outfpath):
                outfpath = untar(taxdump_fpath, self.ncbi_dir, fname)
            logger.debug(f'UPDATE (ncbi,{option}): {outfpath}')
            self.config.set('ncbi', option, outfpath)
        if self.dryrun:
            return
        taxdump_md5 = f'{taxdump_fpath}.md5'
        write_checksum(taxdump_fpath, taxdump_md5)
        os.remove(taxdump_fpath)

    def download_ncbi_files(self, options):
        """Download NCBI database files.

        Parameters
        ----------
        options : iterable
            iterable containing options in 'ncbi' section to download.

        Returns
        -------
        NoneType
            Will update provided `options` in `self.config`.

        Raises
        -------
        ConnectionError
            NCBI file download failed.

        """
        for option in options:
            host = self.config.get('ncbi', 'host')
            ftp_fullpath = self.config.get('database_urls', option)
            ftp_fpath = ftp_fullpath.split(host)[-1]

            if self.config.has_option('ncbi', option) and self.config.get('ncbi', option) is not None:
                outfpath = self.config.get('ncbi', option)
            else:
                outfname = os.path.basename(ftp_fpath)
                outfpath = os.path.join(self.ncbi_dir, outfname)

            logger.debug(f'UPDATE: (ncbi,{option}): {outfpath}')
            self.config.set('ncbi', option, outfpath)

            if self.dryrun:
                return

            with FTP(host) as ftp, open(outfpath, 'wb') as fp:
                ftp.login()
                logger.debug(f'starting {option} download')
                result = ftp.retrbinary(f'RETR {ftp_fpath}', fp.write)
                if not result.startswith('226 Transfer complete'):
                    raise ConnectionError(f'{option} download failed')
                ftp.quit()
        if 'taxdump' in options:
            self.extract_taxdump()
        if 'nr' in options:
            self.format_nr()

    def press_hmms(self):
        """hmmpress markers database files.

        Returns
        -------
        NoneType

        """
        hmms = (os.path.join(self.markers_dir, fn) for fn in os.listdir(self.markers_dir) if fn.endswith('.hmm'))
        hmm_search_str = os.path.join('autometa/databases/markers/','*.h3*')
        pressed_hmms = glob(hmm_search_str)
        pressed_hmms = set(os.path.realpath(os.path.splitext(fp)[0]) for fp in pressed_hmms)
        for fp in hmms:
            if fp in pressed_hmms:
                continue
            hmmer.hmmpress(fp)

    def download_markers(self, options):
        """Download markers database files and amend user config to reflect this.

        Parameters
        ----------
        options : iterable
            iterable containing options in 'markers' section to download.

        Returns
        -------
        NoneType
            Will update provided `options` in `self.config`.

        Raises
        -------
        ConnectionError
            marker file download failed.

        """
        for option in options:
            # First retrieve the markers file url from `option` in `markers` section
            url = self.config.get('database_urls', option)
            if self.config.has_option('markers', option):
                outfpath = self.config.get('markers', option)
            else:
                outfname = os.path.basename(url)
                outfpath = os.path.join(self.markers_dir, outfname)

            if self.dryrun:
                logger.debug(f'UPDATE: (markers,{option}): {outfpath}')
                self.config.set('markers', option, outfpath)
                continue

            # Retrieve markers file and write contents to `outfpath`
            with requests.Session() as session, open(outfpath, 'w') as fh:
                resp = session.get(url)
                if not resp.ok:
                    raise ConnectionError(f'Failed to retrieve {url}')
                fh.write(resp.text)
            self.config.set('markers', option, outfpath)
        self.press_hmms()

    def get_missing(self, section=None, validate=False):
        """Get all missing database files in `options` from `sections`
        in config.

        If the database file already exists, an md5 checksum comparison will be
        performed bewteen the current and file's remote md5 to ensure file integrity
        prior to checking the respective file as valid.

        Parameters
        ----------
        section : str, optional
            Configure provided `section`. Choices include 'markers' and 'ncbi'.
            (default will download/format all database directories)
        validate : bool, optional
            check whether database files are missing and exit (default is False).

        Returns
        -------
        bool or dict
            * if `validate` is True : bool
                all available evaluates to False, otherwise True
            * if `validate` is False : dict
                {section:{option, option,...}, section:{...}, ...}

        """
        sections = [section] if section else Databases.SECTIONS.keys()
        missing = {}
        for section in sections:
            for option in self.config.options(section):
                # Skip user added options not required by Autometa
                if option not in Databases.SECTIONS.get(section):
                    continue
                fpath = self.config.get(section, option)
                if os.path.exists(fpath):
                    continue
                if validate:
                    return True
                if section in missing:
                    missing[section].add(option)
                else:
                    missing.update({section:set([option])})
        # Log missing options
        for section,options in missing.items():
            for option in options:
                logger.debug(f'MISSING: ({section},{option})')
        # ReturnType based on provided `validate` parameter
        if validate:
            return False
        else:
            return missing

    def download_missing(self, section=None):
        """Download missing Autometa database dependencies from provided `section`.
        If no `section` is provided will check all sections.

        Parameters
        ----------
        section : str, optional
            Section to check for missing database files (the default is None).
            Choices include 'ncbi' and 'markers'.

        Returns
        -------
        NoneType
            Will update provided `section` in `self.config`.

        Raises
        -------
        ValueError
            Provided `section` does not match 'ncbi' and 'markers'.

        """
        dispatcher = {'ncbi':self.download_ncbi_files, 'markers':self.download_markers}
        if section and section not in dispatcher:
            raise ValueError(f'{section} does not match "ncbi" or "markers"')
        if section:
            missing = self.get_missing(section=section)
            options = missing.get(section, [])
            dispatcher[section](options)
        else:
            self.get_missing()
            for section,options in missing.items():
                dispatcher[section](options)

    def compare_checksums(self, section=None, validate=False):
        """Get all invalid database files in `options` from `section`
        in config. An md5 checksum comparison will be performed between the
        current and file's remote md5 to ensure file integrity prior to checking
        the respective file as valid.

        Parameters
        ----------
        section : str, optional
            Configure provided `section`. Choices include 'markers' and 'ncbi'.
            (default will download/format all database directories)
        validate : bool, optional
            check whether database files are invalid and exit (default is False).

        Returns
        -------
        bool or dict
            * if `validate` is True : bool
                False if all database files are valid otherwise True.
            * if `validate` is False : dict
                {section:{option, option,...}, section:{...}, ...}

        """
        sections = [section] if section else Databases.SECTIONS.keys()
        invalid = {}
        taxdump_checked = False
        for section in sections:
            for option in self.config.options(section):
                if option not in Databases.SECTIONS.get(section):
                    # Skip user added options not required by Autometa
                    continue
                # nodes.dmp, names.dmp and merged.dmp are all in taxdump.tar.gz
                option = 'taxdump' if option in {'nodes','names','merged'} else option
                fpath = self.config.get(section, option)
                fpath_md5 = f'{fpath}.md5'
                # We can not checksum a file that does not exist.
                if not os.path.exists(fpath) and not os.path.exists(fpath_md5):
                    if validate:
                        return True
                    continue
                # To not waste time checking the taxdump files 3 times.
                if option == 'taxdump' and taxdump_checked:
                    continue
                if os.path.exists(fpath_md5):
                    current_checksum = read_checksum(fpath_md5)
                else:
                    current_checksum = calc_checksum(fpath)
                current_hash, __ = current_checksum.split()
                try:
                    remote_checksum = self.get_remote_checksum(section, option)
                    remote_hash, __ = remote_checksum.split()
                except ConnectionError as err:
                    # Do not mark file as invalid if a connection error occurs.
                    logger.warning(err)
                    continue
                if option == 'taxdump':
                    taxdump_checked = True
                if remote_hash == current_hash:
                    logger.debug(f'{option} checksums match, skipping...')
                    continue
                if section in invalid:
                    invalid[section].add(option)
                else:
                    invalid.update({section:set([option])})
        # Log invalid options
        for section,options in invalid.items():
            for option in options:
                logger.debug(f'INVALID: ({section},{option})')
        # ReturnType based on provided `validate` parameter
        if validate:
            return False
        else:
            return invalid

    def fix_invalid_checksums(self, section=None):
        """Download/Update/Format databases where checksums are out-of-date.

        Parameters
        ----------
        section : str, optional
            Configure provided `section`. Choices include 'markers' and 'ncbi'.
            (default will download/format all database directories)

        Returns
        -------
        NoneType
            Will update provided `options` in `self.config`.

        Raises
        -------
        ConnectionError
            Failed to connect to `section` host site.

        """
        dispatcher = {'ncbi':self.download_ncbi_files, 'markers':self.download_markers}
        if section and section not in dispatcher:
            raise ValueError(f'{section} does not match "ncbi" or "markers"')
        if section:
            invalid = self.compare_checksums(section=section)
            options = invalid.get(section, [])
            dispatcher[section](options)
        else:
            invalid = self.compare_checksums()
            for section,options in invalid.items():
                dispatcher[section](options)

    def configure(self, section=None, no_checksum=False):
        """Configures Autometa's database dependencies by first checking missing
        dependencies then comparing checksums to ensure integrity of files.

        Download and format databases for all options in each section.

        This will only perform the download and formatting if `self.dryrun` is False.
        This will update out-of-date databases if `self.do_update` is True.

        Parameters
        ----------
        section : str, optional
            Configure provided `section`. Choices include 'markers' and 'ncbi'.
            (default will download/format all database directories)
        no_checksum : bool, optional
            Do not perform checksum comparisons (Default is False).

        Returns
        -------
        configparser.ConfigParser
            config with updated options in respective databases sections.

        Raises
        -------
        ValueError
            Provided `section` does not match 'ncbi' or 'markers'.
        ConnectionError
            A connection issue occurred when connecting to NCBI or GitHub.

        """
        self.download_missing(section=section)
        if no_checksum:
            return self.config
        self.fix_invalid_checksums(section=section)
        return self.config

def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format='[%(asctime)s %(levelname)s] %(name)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)

    parser = argparse.ArgumentParser(
        description='Main script to configure Autometa database dependencies.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='By default, with no arguments, will download/format databases into default databases directory.')
    parser.add_argument('--config',
        help='</path/to/input/database.config>', default=DEFAULT_FPATH)
    parser.add_argument('--dryrun',
        help='Log configuration actions but do not perform them.',
        action='store_true', default=False)
    parser.add_argument('--update', help='Update all out-of-date databases.',
        action='store_true', default=False)
    parser.add_argument('--update-markers', help='Update out-of-date markers databases.',
        action='store_true', default=False)
    parser.add_argument('--update-ncbi', help='Update out-of-date ncbi databases.',
        action='store_true', default=False)
    parser.add_argument('--check-dependencies',
        help='Check database dependencies are satisfied.',
        action='store_true', default=False)
    parser.add_argument(
        '--no-checksum',
        help='Do not perform checksum comparisons to validate integrity of database files',
        action='store_true', default=False)
    parser.add_argument('--nproc',
        help='num. cpus to use for DB formatting.', type=int, default=mp.cpu_count())
    parser.add_argument('--out', help='</path/to/output/database.config>')
    args = parser.parse_args()

    config = get_config(args.config)
    dbs = Databases(
        config=config,
        dryrun=args.dryrun,
        nproc=args.nproc,
        do_update=args.update)

    compare_checksums = False
    for update_section in [args.update, args.update_markers, args.update_ncbi]:
        if update_section and not args.no_checksum:
            compare_checksums = True

    if args.update_markers:
        section = 'markers'
    elif args.update_ncbi:
        section = 'ncbi'
    else:
        section = None

    if args.check_dependencies:
        dbs_satisfied = dbs.satisfied(section=section, compare_checksums=compare_checksums)
        logger.info(f'Database dependencies satisfied: {dbs_satisfied}')
        import sys;sys.exit(0)

    # logger.debug(f'Configuring databases')
    config = dbs.configure(section=section, no_checksum=args.no_checksum)

    if not args.out:
        import sys;sys.exit(0)
    put_config(config, args.out)
    logger.info(f'{args.out} written.')

if __name__ == '__main__':
    main()
