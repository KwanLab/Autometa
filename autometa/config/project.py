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

Configuration handling for Autometa User Project.
"""


import logging
import os

from configparser import NoOptionError

from autometa.config import DEFAULT_CONFIG
from autometa.config import get_config
from autometa.config import parse_args
from autometa.config import put_config
from autometa.common.utilities import make_inputs_checkpoints
from autometa.common.utilities import get_existing_checkpoints
from autometa.common.utilities import merge_checkpoints


logger = logging.getLogger(__name__)


class Project:
    """Autometa Project class to configure project directory given `config_fpath`

    Parameters
    ----------
    config_fpath : str
        </path/to/project.config>

    Attributes
    ----------
    dirpath : str
        Path to directory containing `config_fpath`
    config : config.ConfigParser
        interpolated config object parsed from `config_fpath`.
    n_metagenomes : int
        Number of metagenomes contained in project directory
    metagenomes : dict
        metagenomes pertaining to project keyed by number and values of metagenome.config file path.
    new_metagenome_num : int
        Retrieve new minimum metagenome num from metagenomes in project.

    Methods
    ----------
    * self.save()
    * self.new_metagenome_directory()
    * self.setup_checkpoints_and_files()
    * self.add()
    * self.update()
    """

    def __init__(self, config_fpath):
        self.config_fpath = config_fpath
        self.dirpath = os.path.dirname(os.path.realpath(config_fpath))
        self.config = get_config(self.config_fpath)
        if not self.config.has_section("metagenomes"):
            self.config.add_section("metagenomes")

    @property
    def n_metagenomes(self):
        """Return the number of metagenome directories present in the project

        Returns
        -------
        int
            Number of metagenomes contained in project.
        """
        return len(self.metagenomes)

    @property
    def metagenomes(self):
        """retrieve metagenome configs from project.config

        Returns
        -------
        dict
            {metagenome_num:</path/to/metagenome.config>, ...}
        """
        return {
            int(k.strip("metagenome_")): v
            for k, v in self.config.items("metagenomes")
            if os.path.exists(v)
        }

    @property
    def new_metagenome_num(self):
        """Retrieve new minimum metagenome num from metagenomes in project.

        Returns
        -------
        int
            New metagenome number in project.

        """
        # I.e. no metagenomes have been added to project yet.
        if not self.metagenomes:
            return 1
        # max corresponds to highest metagenome number recovered in project directory
        max_num = max(self.metagenomes)
        if max_num == self.n_metagenomes:
            return self.n_metagenomes + 1
        # Otherwise metagenome_num in between max and others has been removed
        # Therefore new metagenome may be inserted.
        for mg_num in range(1, max_num):
            if mg_num in self.metagenomes:
                continue
            return mg_num

    def save(self):
        """Save project config in project directory
        """
        put_config(self.config, self.config_fpath)

    def new_metagenome_directory(self):
        """Create a new metagenome directory in project

        Returns
        -------
        str
            Path to newly created metagenome directory contained in project

        Raises
        ------
        IsADirectoryError
            Directory that is trying to be created already exists
        """
        metagenome_name = f"metagenome_{self.new_metagenome_num:03d}"
        metagenome_dirpath = os.path.join(self.dirpath, metagenome_name)
        # Check presence of metagenome directory
        if os.path.exists(metagenome_dirpath):
            raise IsADirectoryError(metagenome_dirpath)
        os.makedirs(metagenome_dirpath)
        return metagenome_dirpath

    def setup_checkpoints_and_files(self, config, dirpath):
        """Update config files section with symlinks of existing files to metagenome output directory.
        Also get checkpoints from each existing file and write these to a checkpoints file.

        Note
        ----
        Will write checkpoints to `config.get("files", "checkpoints")` file path. Will skip writing checkpoints
        if "checkpoints" is not available in "files".

        Parameters
        ----------
        config : config.ConfigParser
            metagenome config to be updated
        dirpath : str
            Path to output metagenome directory

        Returns
        -------
        config.ConfigParser
            Updated metagenome config
        """
        # symlink any files that already exist and were specified
        checkpoint_inputs = []
        try:
            checkpoints_fpath = config.get("files", "checkpoints")
        except NoOptionError:
            logger.debug("checkpoints option unavailable, skipping.")
            checkpoints_fpath = None
        for option in config.options("files"):
            default_fname = os.path.basename(DEFAULT_CONFIG.get("files", option))
            option_fpath = os.path.realpath(config.get("files", option))
            if os.path.exists(option_fpath):
                if option_fpath.endswith(".gz") and not default_fname.endswith(".gz"):
                    default_fname += ".gz"
                full_fpath = os.path.join(dirpath, default_fname)
                os.symlink(option_fpath, full_fpath)
                checkpoint_inputs.append(full_fpath)
            else:
                full_fpath = os.path.join(dirpath, default_fname)
            config.set("files", option, full_fpath)
        if checkpoints_fpath:
            logger.debug(
                f"Making {len(checkpoint_inputs)} checkpoints and writing to {checkpoints_fpath}"
            )
            checkpoints = make_inputs_checkpoints(checkpoint_inputs)
            checkpoints.to_csv(checkpoints_fpath, sep="\t", index=False, header=True)
        return config

    def add(self, fpath):
        """Setup Autometa metagenome directory given a metagenome.config file.

        Parameters
        ----------
        fpath : str
            </path/to/metagenome.config>

        Returns
        -------
        argparse.Namespace

        Raises
        -------
        IsADirectoryError
            Metagenome output directory already exists
        """
        metagenome_dirpath = self.new_metagenome_directory()
        metagenome_name = os.path.basename(metagenome_dirpath)
        mg_config = get_config(fpath)
        # Add/Update database and env sections for debugging individual metagenome binning runs.
        for section in ["databases", "environ", "versions"]:
            if not mg_config.has_section(section):
                mg_config.add_section(section)
            for option, value in self.config.items(section):
                if not mg_config.has_option(section, option):
                    mg_config.set(section, option, value)
        # symlink any files that already exist and were specified and checkpoint existing files
        self.setup_checkpoints_and_files(config=mg_config, dirpath=metagenome_dirpath)
        # Set outdir parameter and add config section linking metagenome config to project config
        mg_config.set("parameters", "outdir", metagenome_dirpath)
        mg_config_fpath = os.path.join(metagenome_dirpath, f"{metagenome_name}.config")
        mg_config.add_section("config")
        mg_config.set("config", "project", self.config_fpath)
        mg_config.set("config", "metagenome", mg_config_fpath)
        # Save metagenome config to metagenome directory metagenome_00d.config
        put_config(mg_config, mg_config_fpath)
        self.config.set("metagenomes", metagenome_name, mg_config_fpath)
        # Only write updated project config after successful metagenome configuration.
        self.save()
        logger.debug(
            f"updated {self.config_fpath} metagenome: {metagenome_name} : {mg_config_fpath}"
        )
        return parse_args(mg_config_fpath)

    def update(self, metagenome_num, fpath):
        """Update project config metagenomes section with input metagenome.config file.

        Parameters
        ----------
        metagenome_num: int
            metagenome number to update
        fpath : str
            </path/to/new/metagenome.config> This config will overwrite any values in old config
            that are different

        Returns
        -------
        argparse.Namespace

        Raises
        -------
        ValueError
            `metagenome` must be an int and within project config!
        """
        metagenome = f"metagenome_{metagenome_num:03d}"
        if not self.config.has_option("metagenomes", metagenome):
            raise ValueError(
                f"{metagenome_num} must be an int and within project config!"
            )
        old_config_fp = self.config.get("metagenomes", metagenome)
        old_config = get_config(old_config_fp)
        new_config = get_config(fpath)
        new_checkpoints = []
        for section in new_config.sections():
            if not old_config.has_section(section):
                old_config.add_section(section)
            for option in new_config.options(section):
                new_value = new_config.get(section, option)
                if section == "files":
                    if not os.path.exists(new_value):
                        continue
                    new_checkpoints.append(new_value)
                old_config.set(section, option, new_value)
        checkpoints_fpath = old_config.get("files", "checkpoints")
        checkpoints = merge_checkpoints(
            old_checkpoint_fpath=checkpoints_fpath,
            new_checkpoints=new_checkpoints,
            overwrite=True,
        )
        put_config(old_config, old_config_fp)
        logger.debug(f"Updated {metagenome}.config with {fpath}")
        return parse_args(old_config_fp)


def main():
    import argparse
    import logging as logger

    logger.basicConfig(
        format="[%(asctime)s %(levelname)s] %(name)s: %(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
        level=logger.DEBUG,
    )
    parser = argparse.ArgumentParser(
        description="""
    Contains Project class used to manipulate user's Project.
    main logs status of project.
    """
    )
    parser.add_argument("config", help="</path/to/project.config>")
    args = parser.parse_args()
    project = Project(args.config)
    logger.info(
        f"{project.config_fpath} has {project.n_metagenomes} metagenomes in {project.dirpath}"
    )
    logger.info(f"metagenome config numbers: {','.join(map(str,project.metagenomes))}")


if __name__ == "__main__":
    main()
