#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Configuration handling for Autometa User Project.
"""


import copy
import logging
import os
import argparse

from configparser import ConfigParser
from configparser import ExtendedInterpolation

from autometa.config import DEFAULT_FPATH
from autometa.config import DEFAULT_CONFIG
from autometa.config import AUTOMETA_DIR
from autometa.config import get_config
from autometa.config import put_config


logger = logging.getLogger(__name__)


def update_project(config, args):
    """Update input args within config file.

    Parameters
    ----------
    args : argparse.Namespace
        Description of parameter `arg`.

    Returns
    -------
    configparser.ConfigParser
        Updated ConfigParser with user parameters

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    param_section = 'parameters'
    metagenome_section = 'metagenomes'
    config.add_section(param_section)
    for param,value in args.__dict__.items():
        if param == metagenome_section:
            config.add_section(metagenome_section)
            metagenome_count = 1
            for metagenome in value:
                mg = f'metagenome_{metagenome_count:03d}'
                config.set(metagenome_section, mg, metagenome)
                metagenome_count += 1
        if param == 'config' and not value:
            value = DEFAULT_FPATH
        config.set(param_section, param, str(value))
    return config

def setup_project(config):
    """Write project.config to <new project directory>. If directory does not exist,
    one will be made.

    Parameters
    ----------
    config : configparser.ConfigParser
        project config required sections: ['parameters']

    Returns
    -------
    str
        </path/to/projects/project_num/project.config>

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    projects_dp = config.get('parameters','projects')
    try:
        project_num = config.getint('parameters','project')
    except ValueError as err:
        n_projects = 0
        for dp in os.listdir(projects_dp):
            if 'project_' in dp and os.path.isdir(os.path.join(projects_dp,dp)):
                n_projects += 1
        project_num = n_projects + 1
    config.set('parameters','project',str(project_num))
    project_dirname = f'project_{project_num:03d}'
    project_dirpath = os.path.join(projects_dp, project_dirname)
    if not os.path.exists(project_dirpath):
        os.makedirs(project_dirpath)
    config_fp = os.path.join(project_dirpath, 'project.config')
    outconfig = copy.deepcopy(config)
    outconfig.remove_option('parameters','metagenomes')
    outconfig.remove_option('parameters','resume')
    put_config(outconfig, config_fp)
    return config_fp

def setup_metagenome(config):
    """Setup Autometa metagenome directory given a config file of metagenome
    Submission [files] and [parameters] sections.

    Parameters
    ----------
    config : str or configparser.ConfigParser
        </path/to/metagenome.config> or already loaded metagenome ConfigParser

    Returns
    -------
    str
        </path/to/projects/project/metagenome/metagenomeNum.config>

    Raises
    -------
    FileNotFoundError
        project directory or project config does not exist.

    """
    if type(config) is str and os.path.exists(config):
        config = get_config(config)
    mg_config = copy.deepcopy(config)
    # Determine what project metagenome belongs...
    projects_dirpath = mg_config.get('parameters','projects')
    project_num = mg_config.getint('parameters','project')
    project_dname = f"project_{project_num:03d}"
    project_dirpath = os.path.join(projects_dirpath,project_dname)
    if not os.path.exists(project_dirpath):
        raise FileNotFoundError(f'ProjectDirectoryNotFound: {project_dirpath}')
    project_config_fpath = os.path.join(project_dirpath, 'project.config')
    if not os.path.exists(project_config_fpath):
        raise FileNotFoundError(project_config_fpath)
    # Determine metagenome number added to project and update project.config
    metagenomes = [dpath for dpath in os.listdir(project_dirpath)
        if 'metagenome_' in dpath and os.path.isdir(os.path.join(project_dirpath,dpath))]
    metagenome_num = 1 + len(metagenomes)
    metagenome_dirname = f'metagenome_{metagenome_num:03d}'
    metagenome_dirpath = os.path.join(project_dirpath, metagenome_dirname)
    if os.path.exists(metagenome_dirpath):
        raise FileExistsError(metagenome_dirpath)
    os.makedirs(metagenome_dirpath)
    metagenome_fpath = config.get('files','metagenome')
    proj_config = get_config(project_config_fpath)
    proj_config.set('metagenomes', metagenome_dirname, metagenome_fpath)
    put_config(proj_config, project_config_fpath)
    for section in ['databases','environ']:
        if not mg_config.has_section(section):
            mg_config.add_section(section)
        for option,value in proj_config.items(section):
            mg_config.set(section,option,value)
    #Remove project config sections/options if project.config was provided.
    # Change mg config section and parameters to suit respective directory.
    if mg_config.has_section('metagenomes'):
        mg_config.remove_section('metagenomes')
    if mg_config.has_option('parameters','metagenomes'):
        mg_config.remove_option('parameters','metagenomes')
    #symlink any files that already exist and were specified
    for option in mg_config.options('files'):
        default_fname = os.path.basename(DEFAULT_CONFIG.get('files',option))
        fpath = mg_config.get('files',option)
        if os.path.exists(fpath):
            if fpath.endswith('.gz') and not default_fname.endswith('.gz'):
                default_fname += '.gz'
            full_fpath = os.path.join(metagenome_dirpath, default_fname)
            os.symlink(os.path.realpath(fpath),full_fpath)
        elif os.path.basename(fpath).title() == 'None':
            full_fpath = os.path.join(metagenome_dirpath, default_fname)
        else:
            fname = os.path.basename(fpath)
            full_fpath = os.path.join(metagenome_dirpath, fname)
        mg_config.set('files', option, full_fpath)
    mg_config.set('parameters','outdir', metagenome_dirpath)
    mg_config_fpath = os.path.join(metagenome_dirpath, f'{metagenome_dirname}.config')
    mg_config.add_section('config')
    mg_config.set('config','project', project_config_fpath)
    mg_config.set('config','metagenome', mg_config_fpath)
    put_config(mg_config, mg_config_fpath)
    logger.debug(f'updated {project_config_fpath} metagenomes: {metagenome_dirname} : {mg_config_fpath}')
    # Only write updated project config after successful metagenome configuration.
    return mg_config_fpath

def setup_metagenomes(project_config):
    """Build directories for each provided metagenome in `project_config`.

    Parameters
    ----------
    project_config : configparser.ConfigParser
        Description of parameter `project_config`.

    Returns
    -------
    dict
        {
        'metagenome_num':'</path/to/projects/project_num/metagenome_num/metagenome_num.config',
        'metagenome_num':'</path/to/projects/project_num/metagenome_num/metagenome_num.config',
        ...}

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    projects_dirpath = project_config.get('parameters','projects')
    project_dirname = f"project_{project_config.getint('parameters','project'):03d}"
    project_dirpath = os.path.join(projects_dirpath, project_dirname)
    config_fp = os.path.join(project_dirpath, 'project.config')
    config_fpaths = {}
    for metagenome in project_config.options('metagenomes'):
        user_metagenome = project_config.get('metagenomes', metagenome)
        mg_config = copy.deepcopy(project_config)
        mg_config.set('files','metagenome',user_metagenome)
        if mg_config.has_section('metagenomes'):
            mg_config.remove_section('metagenomes')
        written_config = setup_metagenome(config=mg_config)
        config_fpaths.update({metagenome:written_config})
    return config_fpaths

def configure(config, args):
    """Configure output directory with `args`.

    Parameters
    ----------
    args : argparse.Namespace
        Description of parameter `args`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    FileNotFoundError
        config file path does not exist

    """
    config = update_project(config, args)
    return setup_project(config)

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
