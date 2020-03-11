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

AutometaUser configuration class
"""


import logging
import os

import argparse

# TODO: Refactor autometa.config later as AutometaConfigUtils lib or something
from autometa.config import get_config
from autometa.config import put_config
from autometa.config import parse_config
from autometa.config import AUTOMETA_DIR
from autometa.config import DEFAULT_CONFIG
from autometa.config import DEFAULT_FPATH
from autometa.config import databases
from autometa.config import environ
from autometa.config.project import Project
from autometa.common import utilities


logger = logging.getLogger(__name__)


class AutometaUser:
    """docstring for AutometaUser."""

    def __init__(self, config_fpath=None, dryrun=True, nproc=2):
        self.dryrun= dryrun
        self.nproc = nproc
        self.config_fp = config_fpath
        self.config = get_config(self.config_fp) if self.config_fp else DEFAULT_CONFIG
        if not self.config.has_section('common'):
            self.config.add_section('common')
        self.config.set('common','home_dir', AUTOMETA_DIR)

    def configure(self, configure_environ=True, configure_databases=True):
        if configure_environ:
            self.config = environ.configure(self.config)
        if configure_databases:
            self.config = databases.configure(self.config, dryrun=self.dryrun, nproc=self.nproc)

    def new_workspace(self, fpath):
        """Configure new project at `outdir`.

        Parameters
        ----------
        fpath : str
            </path/to/workspace/project_<num>/project.config>

        Returns
        -------
        autometa.config.project.Project object

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        # 1. configure project from default config and provided config file
        self.configure()
        dpath = os.path.dirname(fpath)
        if not os.path.exists(dpath):
            os.makedirs(dpath)
        put_config(self.config, fpath)
        return Project(fpath)

    def prepare_run(self, config_fpath):
        """Prepares metagenome binning run using provided `config_fpath`.

        This method performs a number of configuration checks to ensure the
        binning run will perform without conflicts.
        1. workspace check: Will construct workspace directory if provided does not
        exist.
        2. Project check: Will configure a new project if project number is not
        found in workspace directory.
        3. Metagenome check: Will update if existing with edits or resume
        if existing without edits. Otherwise will add new metagenome to project.

        Parameters
        ----------
        config_fpath : str
            </path/to/metagenome.config>

        Returns
        -------
        argparse.Namespace
            access to parameters and files from config via syntax...
            i.e.
            generate namespace:
                mgargs = prepare_run(mg_config)
            access namespace:
                mgargs.files.<file>
                mgargs.parameters.<parameter>

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        mgargs = parse_config(config_fpath)
        # 1 check workspace exists
        workspace = os.path.realpath(mgargs.parameters.workspace)
        if not os.path.exists(workspace):
            os.makedirs(workspace)
        # 2 check project exists
        proj_name = f'project_{mgargs.parameters.project:03d}'
        project_dirpath = os.path.realpath(os.path.join(workspace,proj_name))
        project_config_fp = os.path.join(project_dirpath, 'project.config')
        if not os.path.exists(project_dirpath) or not os.path.exists(project_config_fp):
            project = self.new_workspace(project_config_fp)
        else:
            project = Project(project_config_fp)
        # 3 check whether existing or new run with metagenome_num
        metagenome = f'metagenome_{mgargs.parameters.metagenome_num:03d}'
        if metagenome not in project.metagenomes:
            mgargs = project.add(config_fpath)
            project.save()
            return mgargs
        # If resuming existing metagenome run. Check whether config file has changed.
        old_config_fp = project.metagenomes.get(metagenome)
        old_chksum = utilities.get_checksum(old_config_fp)
        new_chksum = utilities.get_checksum(config_fpath)
        if old_chksum != new_chksum:
            mgargs = project.update(
                metagenome_num=mgargs.parameters.metagenome_num,
                fpath=config_fpath)
        project.save()
        return mgargs

def main(args):
    logger.info(args.user)

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Concise Functional Description of Script')
    parser.add_argument('user',help='</path/to/user.config>')
    args = parser.parse_args()
    main(args)
