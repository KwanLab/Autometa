#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
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

Configuration handling for Autometa User Project.
"""


import copy
import logging
import os
import argparse

from autometa.config import DEFAULT_CONFIG
from autometa.config import get_config
from autometa.config import parse_config
from autometa.config import put_config


logger = logging.getLogger(__name__)


class Project:
    """Autometa Project."""

    def __init__(self, config_fpath):
        self.config_fpath = config_fpath
        self.dirpath = os.path.dirname(os.path.realpath(config_fpath))
        self.config = get_config(self.config_fpath)
        if not self.config.has_section('metagenomes'):
            self.config.add_section('metagenomes')

    @property
    def n_metagenomes(self):
        return len(self.metagenomes)

    @property
    def metagenomes(self):
        """retrieve metagenome configs from project.config

        Returns
        -------
        dict
            {metagenome_num:</path/to/metagenome.config>, ...}
        """
        return {int(k.strip('metagenome_')):v for k,v in self.config.items('metagenomes') if os.path.exists(v)}

    def new_metagenome_num(self):
        """Retrieve new minimum metagenome num from metagenomes in project.

        Returns
        -------
        int
            Description of returned object.

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        # I.e. no metagenomes have been added to project yet.
        if not self.metagenomes:
            return 1
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
        put_config(self.config, self.config_fpath)

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
        FileNotFoundError
            Directory found but metagenome.config not present
        IsADirectoryError
            Metagenome output directory already exists
        """
        # metagenome_num = 1 + self.n_metagenomes
        metagenome_num = self.new_metagenome_num()
        metagenome_name = f'metagenome_{metagenome_num:03d}'
        metagenome_dirpath = os.path.join(self.dirpath, metagenome_name)
        mg_config_fpath = os.path.join(metagenome_dirpath, f'{metagenome_name}.config')
        # Check presence of metagenome directory and config
        mg_config_present = os.path.exists(mg_config_fpath)
        mg_dir_present = os.path.exists(metagenome_dirpath)
        if not mg_config_present and mg_dir_present:
            raise FileNotFoundError(f'{mg_config_fpath} is not present but the directory exists! Either remove the directory or locate the config file before continuing.')
        if mg_dir_present:
            raise IsADirectoryError(metagenome_dirpath)

        os.makedirs(metagenome_dirpath)
        mg_config = get_config(fpath)
        # Add database and env for debugging individual metagenome binning runs.
        for section in ['databases','environ','versions']:
            if not mg_config.has_section(section):
                mg_config.add_section(section)
            for option,value in self.config.items(section):
                mg_config.set(section,option,value)
        #symlink any files that already exist and were specified
        for option in mg_config.options('files'):
            default_fname = os.path.basename(DEFAULT_CONFIG.get('files',option))
            option_fpath = os.path.realpath(mg_config.get('files',option))
            if os.path.exists(option_fpath):
                if option_fpath.endswith('.gz') and not default_fname.endswith('.gz'):
                    default_fname += '.gz'
                full_fpath = os.path.join(metagenome_dirpath, default_fname)
                os.symlink(option_fpath,full_fpath)
            else:
                full_fpath = os.path.join(metagenome_dirpath, default_fname)
            mg_config.set('files', option, full_fpath)
        mg_config.set('parameters','outdir', metagenome_dirpath)
        mg_config_fpath = os.path.join(metagenome_dirpath, f'{metagenome_name}.config')
        mg_config.add_section('config')
        mg_config.set('config','project', self.config_fpath)
        mg_config.set('config','metagenome', mg_config_fpath)
        put_config(mg_config, mg_config_fpath)
        # Only write updated project config after successful metagenome configuration.
        self.config.set('metagenomes',metagenome_name,mg_config_fpath)
        logger.debug(f'updated {self.config_fpath} metagenome: {metagenome_name} : {mg_config_fpath}')
        return parse_config(mg_config_fpath)

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
        metagenome = f'metagenome_{metagenome_num:03d}'
        if not self.config.has_option('metagenomes',metagenome):
            raise ValueError(f'{metagenome_num} must be an int and within project config!')
        old_config_fp = self.config.get('metagenomes', metagenome)
        old_config = get_config(old_config_fp)
        new_config = get_config(fpath)
        for section in new_config.sections():
            if not old_config.has_section(section):
                old_config.add_section(section)
            for option in new_config.options(section):
                new_value = new_config.get(section,option)
                # TODO: Update file checksums (checkpoint update)
                # TODO: Check if new value exists... Otherwise keep old option
                if section == 'files' and not os.path.exists(new_value):
                    continue
                old_config.set(section, option, new_value)
        put_config(old_config, old_config_fp)
        logger.debug(f'Updated {metagenome}.config with {fpath}')
        return parse_config(old_config_fp)

if __name__ == '__main__':
    import sys;sys.exit(1)
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Autometa Project configuration')
    parser.add_argument('--config',help='<help text of positional arg>')
    args = parser.parse_args()
    main(args)
