#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Autometa User Configuration Class
"""


import logging
import os

import argparse

from autometa.config import get_config
from autometa.config import parse_config
from autometa.config import DEFAULT_CONFIG
from autometa.config import AUTOMETA_DIR
from autometa.config import databases
from autometa.config import environ
from autometa.config import project


logger = logging.getLogger(__name__)


class AutometaUser:
    """docstring for AutometaUser."""

    def __init__(self, config_fpath=None, dryrun=True):
        self.dryrun= dryrun
        self.config_fp = config_fpath
        self.config = get_config(self.config_fp) if self.config_fp else DEFAULT_CONFIG
        if not self.config.has_section('common'):
            self.config.add_section('common')
        self.config.set('common','home_dir', AUTOMETA_DIR)

    def configure(self, configure_environ=True, configure_databases=True, nproc=2):
        if configure_environ:
            self.config = environ.configure(self.config)
        if configure_databases:
            self.config = databases.configure(self.config, dryrun=self.dryrun, nproc=nproc)

    def new_project(self, args):
        """Configure new project with input args.

        Parameters
        ----------
        args : argparse.Namespace
            Description of parameter `args`.

        Returns
        -------
        dict
            {'project':</path/to/projects/project_num/project.config>,
            'metagenomes':{
                    'metagenome_num':'</path/to/projects/project_num/metagenome_num/metagenome_num.config>',
                    'metagenome_num':'</path/to/projects/project_num/metagenome_num/metagenome_num.config>',
                    ...
                },
            }

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        proj_config = project.configure(self.config, args)
        metagenomes_configs = project.setup_metagenomes(get_config(proj_config))
        return {'project':proj_config, 'metagenomes':metagenomes_configs}

    def add_metagenomes(self, metagenomes_configs):
        mg_configs = {}
        for metagenome_config in metagenomes_configs:
            mg_config_fpath = project.setup_metagenome(metagenome_config)
            mg_name = os.path.basename(mg_config_fpath).strip('.config')
            mgargs = parse_config(mg_config_fpath)
            mg_configs.update({mg_name:mgargs})
        return mg_configs

    def get_mgargs(self, projects_dir, project_num, metagenome_num):
        for arg in [project_num, metagenome_num]:
            if type(arg) is not int:
                raise TypeError(f'{args} is type: {type(arg)}')
        if project_num <= 0:
            raise ValueError(f'project num: {project_num} is invalid')
        if metagenome_num <= 0:
            raise ValueError(f'metagenome_num {metagenome_num} is invalid')
        project_name = f'project_{project_num:03d}'
        metagenome_dirname = f'metagenome_{metagenome_num:03d}'
        metagenome_fname = f'{metagenome_dirname}.config'
        metagenome_config = os.path.join(projects_dir, project_name, metagenome_dirname, metagenome_fname)
        return {metagenome_dirname:parse_config(metagenome_config)}

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
