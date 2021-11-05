#!/usr/bin/env python
"""
COPYRIGHT
Copyright 2021 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
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
import socket
import subprocess
import tempfile

import multiprocessing as mp

from typing import Dict, Iterable
from configparser import ConfigParser
from ftplib import FTP
from glob import glob

from autometa.common.utilities import untar
from autometa.common.utilities import calc_checksum
from autometa.common.utilities import read_checksum
from autometa.common.utilities import write_checksum
from autometa.common.utilities import internet_is_connected
from autometa.common.exceptions import ChecksumMismatchError
from autometa.common.external import diamond
from autometa.common.external import hmmscan
from autometa.config.utilities import DEFAULT_FPATH
from autometa.config.utilities import DEFAULT_CONFIG
from autometa.config.utilities import AUTOMETA_DIR
from autometa.config.utilities import put_config, get_config


logger = logging.getLogger(__name__)
urllib_logger = logging.getLogger("urllib3")
urllib_logger.setLevel(logging.WARNING)


class Databases:
    """Database class containing methods to allow downloading/formatting/updating
    Autometa database dependencies.

    Parameters
    ----------
    config : config.ConfigParser
        Config containing database dependency information.
        (the default is DEFAULT_CONFIG).

    dryrun : bool
        Run through database checking without performing
        downloads/formatting (the default is False).

    nproc : int
        Number of processors to use to perform database formatting.
        (the default is mp.cpu_count()).

    update : bool
        Overwrite existing databases with more up-to-date database
        files. (the default is False).

    Attributes
    ----------
    ncbi_dir : str </path/to/databases/ncbi> markers_dir : str
        </path/to/databases/markers> SECTIONS : dict keys are `sections`
        respective to database config sections and values are options within
        the `sections`.
    """

    SECTIONS = {
        "ncbi": ["nodes", "names", "merged", "accession2taxid", "nr"],
        "markers": [
            "bacteria_single_copy",
            "bacteria_single_copy_cutoffs",
            "archaea_single_copy",
            "archaea_single_copy_cutoffs",
        ],
    }

    def __init__(
        self,
        config=DEFAULT_CONFIG,
        dryrun=False,
        nproc=mp.cpu_count(),
        update=False,
    ):
        """

        At instantiation of Databases instance, if any of the respective
        database directories do not exist, they will be created. This will be
        reflected in the provided `config`.

        Parameters
        ----------
        config : config.ConfigParser
            Config containing database dependency information.
            (the default is DEFAULT_CONFIG).
        dryrun : bool
            Run through database checking without performing downloads/
            formatting (the default is False).

        nproc : int
            Number of processors to use to perform database formatting.
            (the default is mp.cpu_count()).
        update : bool
            Overwrite existing databases with more up-to-date database files.
            (the default is False).

        Returns
        -------
        databases.Databases: Databases object
            instance of Databases class

        """
        if not isinstance(config, ConfigParser):
            raise TypeError(f"config is not ConfigParser : {type(config)}")
        if not isinstance(dryrun, bool):
            raise TypeError(f"dryrun must be boolean. type: {type(dryrun)}")

        self.config = config
        self.dryrun = dryrun
        self.nproc = nproc
        self.update = update
        if self.config.get("common", "home_dir") == "None":
            # neccessary if user is running autometa-update-database
            # before running the autometa-config endpoint.
            # where :func:`~autometa.config.utilities.set_home_dir` would've been called.
            self.config.set("common", "home_dir", AUTOMETA_DIR)
        if not self.config.has_section("databases"):
            self.config.add_section("databases")
        for section in Databases.SECTIONS:
            if self.config.has_option("databases", section):
                continue
            outdir = DEFAULT_CONFIG.get("databases", section)
            self.config.set("databases", section, outdir)
        self.ncbi_dir = self.config.get("databases", "ncbi")
        self.markers_dir = self.config.get("databases", "markers")
        for outdir in {self.ncbi_dir, self.markers_dir}:
            if not os.path.exists(outdir):
                os.makedirs(outdir)

    def satisfied(self, section: str = None, compare_checksums: bool = False) -> bool:
        """Determines whether all database dependencies are satisfied.

        Parameters
        ----------
        section : str
            section to retrieve for `checksums` section.
            Choices include: 'ncbi' and 'markers'.
        compare_checksums : bool, optional
            Also check if database information is up-to-date with current
            hosted databases. (default is False).

        Returns
        -------
        bool
            True if all database dependencies are satisfied, otherwise False.

        """
        any_missing = self.get_missing(section=section)
        if compare_checksums:
            any_invalid = self.compare_checksums(section=section)
        else:
            any_invalid = {}
        return not any_missing and not any_invalid

    def get_remote_checksum(self, section: str, option: str) -> str:
        """Get the checksum from provided `section` respective to `option` in
        `self.config`.

        Parameters
        ----------
        section : str
            section to retrieve for `checksums` section.
            Choices include: 'ncbi' and 'markers'.
        option : str
            `option` in `checksums` section corresponding to the section
            checksum file.

        Returns
        -------
        str
            checksum of remote md5 file. e.g. 'hash filename\n'

        Raises
        -------
        ValueError
            'section' must be 'ncbi' or 'markers'
        ConnectionError
            No internet connection available.
        ConnectionError
            Failed to connect to host for provided `option`.

        """
        if section not in {"ncbi", "markers"}:
            raise ValueError(
                f"'section' must be 'ncbi' or 'markers'. Provided: {section}"
            )
        if not internet_is_connected():
            raise ConnectionError("Cannot connect to the internet")
        if section == "ncbi":
            host = self.config.get(section, "host")
            ftp_fullpath = self.config.get("checksums", option)
            chksum_fpath = ftp_fullpath.split(host)[-1]
            with FTP(host) as ftp, tempfile.TemporaryFile() as fp:
                ftp.login()
                result = ftp.retrbinary(f"RETR {chksum_fpath}", fp.write)
                if not result.startswith("226 Transfer complete"):
                    raise ConnectionError(f"{chksum_fpath} download failed")
                ftp.quit()
                fp.seek(0)
                checksum = fp.read().decode()
        elif section == "markers":
            url = self.config.get("checksums", option)
            with requests.Session() as session:
                resp = session.get(url)
                if not resp.ok:
                    raise ConnectionError(f"Failed to retrieve {url}")
                checksum = resp.text
        return checksum

    def format_nr(self) -> None:
        """Construct a diamond formatted database (nr.dmnd) from `nr` option
        in `ncbi` section in user config.

        NOTE: The checksum 'nr.dmnd.md5' will only be generated if nr.dmnd
        construction is successful. If the provided `nr` option in `ncbi` is
        'nr.gz' the database will be removed after successful database
        formatting.

        Returns
        -------
        NoneType
            config updated option:'nr' in section:'ncbi'.

        """
        db_infpath = self.config.get("ncbi", "nr")
        db_infpath_md5 = f"{db_infpath}.md5"
        db_outfpath = db_infpath.replace(".gz", ".dmnd")

        db_outfpath_exists = os.path.exists(db_outfpath)
        if db_outfpath_exists:
            db_outfpath_hash, __ = calc_checksum(db_outfpath).split()

        remote_checksum_matches = False
        current_nr_checksum_matches = False
        # Check database and database checksum is up-to-date
        if os.path.exists(db_infpath_md5) and db_outfpath_exists:
            # Check if the current db md5 is up-to-date with the remote db md5
            current_hash, __ = read_checksum(db_infpath_md5).split()
            remote_hash, __ = self.get_remote_checksum("ncbi", "nr").split()
            if remote_hash == current_hash:
                remote_checksum_matches = True
            # Check if the current db md5 matches the calc'd db checksum
            if db_outfpath_hash == current_hash:
                current_nr_checksum_matches = True

        db_outfpath_md5 = f"{db_outfpath}.md5"
        db_outfpath_md5_checksum_matches = False
        if os.path.exists(db_outfpath_md5) and db_outfpath_exists:
            db_outfpath_md5_hash, __ = read_checksum(db_outfpath_md5).split()
            if db_outfpath_hash == db_outfpath_md5_hash:
                db_outfpath_md5_checksum_matches = True

        checksum_checks = ["nr.dmnd.md5", "nr.gz.md5", "remote nr.gz.md5"]
        checksum_matches = [
            db_outfpath_md5_checksum_matches,
            current_nr_checksum_matches,
            remote_checksum_matches,
        ]
        for checksum_match, checksum_check in zip(checksum_matches, checksum_checks):
            # If the checksums do not match, we need to update the database file.
            if checksum_match:
                logger.debug(f"{checksum_check} checksum matches, skipping...")
                self.config.set("ncbi", "nr", db_outfpath)
                logger.debug(f"set ncbi nr: {db_outfpath}")
                return
            # Only update out-of-date db files if user wants to update via self.update
            if not self.update and checksum_check == "remote nr.gz.md5":
                return

        diamond.makedatabase(fasta=db_infpath, database=db_outfpath, nproc=self.nproc)
        # Write checksum for nr.dmnd
        write_checksum(db_outfpath, db_outfpath_md5)

        if os.path.basename(db_infpath) == "nr.gz":
            # nr.gz will be removed after successful nr.dmnd construction
            os.remove(db_infpath)

        self.config.set("ncbi", "nr", db_outfpath)
        logger.debug(f"set ncbi nr: {db_outfpath}")

    def extract_taxdump(self) -> None:
        """Extract autometa required files from ncbi taxdump.tar.gz archive
        into ncbi databases directory and update user config with extracted
        paths.

        This only extracts nodes.dmp, names.dmp and merged.dmp from
        taxdump.tar.gz if the files do not already exist. If `update`
        was originally supplied as `True` to the Databases instance, then the
        previous files will be replaced by the new taxdump files.

        After successful extraction of the files, a checksum will be written
        of the archive for future checking.

        Returns
        -------
        NoneType
            Will update `self.config` section `ncbi` with options 'nodes',
            'names','merged'

        """
        taxdump_fpath = self.config.get("ncbi", "taxdump")
        taxdump_files = [
            ("nodes", "nodes.dmp"),
            ("names", "names.dmp"),
            ("merged", "merged.dmp"),
        ]
        for option, fname in taxdump_files:
            outfpath = os.path.join(self.ncbi_dir, fname)
            if self.dryrun:
                logger.debug(f"UPDATE (ncbi,{option}): {outfpath}")
                self.config.set("ncbi", option, outfpath)
                continue
            # Only update the taxdump files if the user says to do an update.
            if self.update and os.path.exists(outfpath):
                os.remove(outfpath)
            # Only extract the taxdump files if this is not a "dryrun"
            if not os.path.exists(outfpath):
                outfpath = untar(taxdump_fpath, self.ncbi_dir, fname)
            write_checksum(outfpath, f"{outfpath}.md5")

            logger.debug(f"UPDATE (ncbi,{option}): {outfpath}")
            self.config.set("ncbi", option, outfpath)

    def download_ncbi_files(self, options: Iterable) -> None:
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
        subprocess.CalledProcessError
            NCBI file download with rsync failed.
        ConnectionError
            NCBI file checksums do not match after file transfer.

        """
        # s.t. set methods are available
        options = set(options)
        # If any of the taxdump.tar.gz files are missing,
        # we need to check that taxdump tarball is available to extract them (see self.extract_taxdump).
        for taxdump_option in {"nodes", "names", "merged"}:
            if taxdump_option in options:
                options.add("taxdump")
                options.discard(taxdump_option)
        for option in options:
            ftp_fullpath = self.config.get("database_urls", option)

            if (
                self.config.has_option("ncbi", option)
                and self.config.get("ncbi", option) is not None
            ):
                outfpath = self.config.get("ncbi", option)
            else:
                outfname = os.path.basename(ftp_fullpath)
                outfpath = os.path.join(self.ncbi_dir, outfname)

            logger.debug(f"UPDATE: (ncbi,{option}): {outfpath}")
            self.config.set("ncbi", option, outfpath)

            if self.dryrun:
                return

            rsync_fpath = ftp_fullpath.replace("ftp", "rsync", 1)
            cmd = ["rsync", "--quiet", "--archive", rsync_fpath, outfpath]
            logger.debug(f"starting {option} download")
            subprocess.run(
                cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
            )
            checksum_outfpath = f"{outfpath}.md5"
            write_checksum(outfpath, checksum_outfpath)
            current_checksum = read_checksum(checksum_outfpath)
            current_hash, __ = current_checksum.split()
            remote_checksum = self.get_remote_checksum("ncbi", option)
            remote_hash, __ = remote_checksum.split()
            if current_hash != remote_hash:
                raise ChecksumMismatchError(f"{option} download failed")
        if "taxdump" in options:
            self.extract_taxdump()
        if "nr" in options:
            self.format_nr()

    def press_hmms(self) -> None:
        """hmmpress markers hmm database files.

        Returns
        -------
        NoneType

        """
        hmm_search_str = os.path.join(self.markers_dir, "*.h3?")
        # First search for pressed hmms to remove from list to hmmpress
        pressed_hmms = {
            os.path.realpath(os.path.splitext(fp)[0])
            for fp in glob(hmm_search_str)
            if not fp.endswith(".md5")
        }
        # Now retrieve all hmms in markers directory
        hmms = (
            os.path.join(self.markers_dir, fn)
            for fn in os.listdir(self.markers_dir)
            if fn.endswith(".hmm")
        )
        # Filter by hmms not already pressed
        hmms = (fpath for fpath in hmms if fpath not in pressed_hmms)
        # Press hmms and write checksums of their indices
        for hmm_fp in hmms:
            hmmscan.hmmpress(hmm_fp)
            for index_fp in glob(f"{hmm_fp}.h3?"):
                write_checksum(index_fp, f"{index_fp}.md5")

    def download_markers(self, options: Iterable) -> None:
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
            # First retrieve the markers file url from `option` in `markers`
            url = self.config.get("database_urls", option)
            if self.config.has_option("markers", option):
                outfpath = self.config.get("markers", option)
            else:
                outfname = os.path.basename(url)
                outfpath = os.path.join(self.markers_dir, outfname)

            if self.dryrun:
                logger.debug(f"UPDATE: (markers,{option}): {outfpath}")
                self.config.set("markers", option, outfpath)
                continue

            # Retrieve markers file and write contents to `outfpath`
            with requests.Session() as session, open(outfpath, "w") as fh:
                resp = session.get(url)
                if not resp.ok:
                    raise ConnectionError(f"Failed to retrieve {url}")
                fh.write(resp.text)
            self.config.set("markers", option, outfpath)
            checksum_outfpath = f"{outfpath}.md5"
            write_checksum(outfpath, checksum_outfpath)
            current_checksum = read_checksum(checksum_outfpath)
            current_hash, __ = current_checksum.split()
            remote_checksum = self.get_remote_checksum("markers", option)
            remote_hash, __ = remote_checksum.split()
            if current_hash != remote_hash:
                raise ChecksumMismatchError(f"{option} download failed")
        self.press_hmms()

    def get_missing(self, section: str = None) -> Dict[str, Dict]:
        """Get all missing database files in `options` from `sections`
        in config.

        Parameters
        ----------
        section : str, optional
            Configure provided `section`. Choices include 'markers' and 'ncbi'.
            (default will download/format all database directories)

        Returns
        -------
        dict
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
                if section in missing:
                    missing[section].add(option)
                else:
                    missing.update({section: set([option])})
        # Log missing options
        for section, options in missing.items():
            for option in options:
                logger.debug(f"MISSING: ({section},{option})")
        return missing

    def download_missing(self, section: str = None) -> None:
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
        dispatcher = {
            "ncbi": self.download_ncbi_files,
            "markers": self.download_markers,
        }
        if section and section not in dispatcher:
            raise ValueError(f"{section} does not match {dispatcher.keys()}")
        if section:
            missing = self.get_missing(section=section)
            options = missing.get(section, [])
            dispatcher[section](options)
        else:
            missing = self.get_missing()
            for section, options in missing.items():
                dispatcher[section](options)

    def compare_checksums(self, section: str = None) -> Dict[str, Dict]:
        """Get all invalid database files in `options` from `section`
        in config. An md5 checksum comparison will be performed between the
        current and file's remote md5 to ensure file integrity prior to
        checking the respective file as valid.

        Parameters
        ----------
        section : str, optional Configure provided `section` Choices include
            'markers' and 'ncbi'. (default will download/format all database
            directories)

        Returns
        -------
        dict {section:{option, option,...}, section:{...}, ...}

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
                option = "taxdump" if option in {"nodes", "names", "merged"} else option
                fpath = self.config.get(section, option)
                fpath_md5 = f"{fpath}.md5"
                # We can not checksum a file that does not exist.
                if not os.path.exists(fpath) and not os.path.exists(fpath_md5):
                    continue
                # To not waste time checking the taxdump files 3 times.
                if option == "taxdump" and taxdump_checked:
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
                if option == "taxdump":
                    taxdump_checked = True
                if remote_hash == current_hash:
                    logger.debug(f"{option} checksums match, skipping...")
                    continue
                if section in invalid:
                    invalid[section].add(option)
                else:
                    invalid.update({section: set([option])})
        # Log invalid options
        for section, options in invalid.items():
            for option in options:
                logger.debug(f"INVALID: ({section},{option})")
        return invalid

    def fix_invalid_checksums(self, section: str = None) -> None:
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
        dispatcher = {
            "ncbi": self.download_ncbi_files,
            "markers": self.download_markers,
        }
        if section and section not in dispatcher:
            raise ValueError(f'{section} does not match "ncbi" or "markers"')
        if section:
            invalid = self.compare_checksums(section=section)
            options = invalid.get(section, set())
            dispatcher[section](options)
        else:
            invalid = self.compare_checksums()
            for section, options in invalid.items():
                dispatcher[section](options)

    def configure(self, section: str = None, no_checksum: bool = False) -> ConfigParser:
        """Configures Autometa's database dependencies by first checking missing
        dependencies then comparing checksums to ensure integrity of files.

        Download and format databases for all options in each section.

        This will only perform the download and formatting if `self.dryrun` is
        False. This will update out-of-date databases if `self.update` is
        True.

        Parameters
        ----------
        section : str, optional Configure provided `section`. Choices include
            'markers' and 'ncbi'. (default will download/format all database
            directories) no_checksum : bool, optional Do not perform checksum
            comparisons (Default is False).

        Returns
        -------
        configparser.ConfigParser config with updated options in respective
            databases sections.

        Raises
        -------
        ValueError Provided `section` does not match 'ncbi' or 'markers'.
            ConnectionError A connection issue occurred when connecting to NCBI
            or GitHub.

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
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )

    parser = argparse.ArgumentParser(
        description="Main script to configure Autometa database dependencies.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="By default, with no arguments, will download/format databases "
        "into default databases directory.",
    )
    parser.add_argument(
        "--config",
        help="</path/to/input/database.config>",
        default=DEFAULT_FPATH,
    )
    parser.add_argument(
        "--dryrun",
        help="Log configuration actions but do not perform them.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--update-all",
        help="Update all out-of-date databases.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--update-markers",
        help="Update out-of-date markers databases.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--update-ncbi",
        help="Update out-of-date ncbi databases.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--check-dependencies",
        help="Check database dependencies are satisfied.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--no-checksum",
        help="Do not perform remote checksum comparisons to validate databases"
        " are up-to-date.",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--nproc",
        help="num. cpus to use for DB formatting.",
        type=int,
        default=mp.cpu_count(),
    )
    parser.add_argument("--out", help="</path/to/output/database.config>")
    args = parser.parse_args()

    config = get_config(args.config)
    dbs = Databases(
        config=config,
        dryrun=args.dryrun,
        nproc=args.nproc,
        update=args.update_all,
    )

    compare_checksums = False
    for update_section in [args.update_all, args.update_markers, args.update_ncbi]:
        # If any of these flags are provided and the --no-checksum flag was NOT provided.
        # then we will go through checksum comparisons for all sections
        if update_section and not args.no_checksum:
            compare_checksums = True

    if args.update_markers:
        section = "markers"
    elif args.update_ncbi:
        section = "ncbi"
    else:
        section = None

    if args.check_dependencies:
        dbs_satisfied = dbs.satisfied(
            section=section, compare_checksums=compare_checksums
        )
        logger.info(f"Database dependencies satisfied: {dbs_satisfied}")
        import sys

        sys.exit(0)

    config = dbs.configure(section=section, no_checksum=args.no_checksum)

    if not args.out:
        import sys

        sys.exit(0)
    put_config(config, args.out)
    logger.info(f"{args.out} written.")


if __name__ == "__main__":
    main()
