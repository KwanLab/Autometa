#!/usr/bin/env python3
"""
Configuration handling for Autometa Databases.
"""


import logging
import os
import requests

from configparser import ConfigParser
from configparser import ExtendedInterpolation
from ftplib import FTP

from autometa.config import get_config
from autometa.config import DEFAULT_CONFIG
from autometa.config import put_config
from autometa.config import AUTOMETA_DIR
from autometa.common.utilities import untar
from autometa.common.external import diamond


logger = logging.getLogger(__name__)

DB_SECTIONS = {
    'ncbi':[
        'nodes',
        'names',
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


def format_nr(config, dryrun, nproc=2):
    outdir = config.get('databases','ncbi')
    nr = config.get('ncbi','nr')
    formatted_nr = os.path.splitext(os.path.basename(nr))[0]
    outfpath = os.path.join(outdir, formatted_nr)
    if not dryrun:
        diamond.makedatabase(fasta=nr, database=outfpath, nproc=nrpoc)
    config.set('ncbi','nr',outfpath)
    logger.debug(f'set ncbi nr to {outfpath}')
    return config

def extract_taxdump(config, dryrun):
    """Short summary.

    Parameters
    ----------
    config : type
        Description of parameter `config`.
    dryrun : type
        Description of parameter `dryrun`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    outdir = config.get('databases','ncbi')
    taxdump = config.get('ncbi','taxdump')
    extraction_files = [('nodes','nodes.dmp'),('names','names.dmp'),('merged','merged.dmp')]
    for option,fname in extraction_files:
        outfp = os.path.join(outdir,fname)
        if not dryrun:
            outfp = untar(taxdump, outdir, fname)
        logger.debug(f'update ncbi : {option} : {outfp}')
        config.set('ncbi',option,outfp)
    return config

def update_ncbi(config, options, dryrun, nproc=2):
    """Update NCBI database.

    Parameters
    ----------
    config : type
        Description of parameter `config`.
    options : type
        Description of parameter `options`.
    dryrun : type
        Description of parameter `dryrun`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    section = 'ncbi'
    if not config.has_section('databases'):
        config.add_section('databases')
    if not config.has_option('databases',section):
        outdir = DEFAULT_CONFIG.get('databases', section)
        config.set('databases', section, outdir)
    else:
        outdir = config.get('databases',section)
    host = DEFAULT_CONFIG.get(section,'host')
    if dryrun:
        for option in options:
            ftp_fullpath = DEFAULT_CONFIG.get('database_urls',option)
            ftp_fpath = ftp_fullpath.split(host)[-1]
            if config.has_option(section, option):
                outfpath = config.get(section, option)
            else:
                outfname = os.path.basename(ftp_fpath)
                outfpath = os.path.join(outdir, outfname)
            logger.debug(f'update {section} : {option} : {outfpath}')
            config.set(section, option, outfpath)
            continue
        return extract_taxdump(config, dryrun)
    with FTP(host) as ftp:
        ftp.login()
        for option in options:
            ftp_fullpath = DEFAULT_CONFIG.get('database_urls',option)
            ftp_fpath = ftp_fullpath.split(host)[-1]
            if config.has_option(section, option):
                outfpath = config.get(section, option)
            else:
                outfname = os.path.basename(ftp_fpath)
                outfpath = os.path.join(outdir, outfname)
            logger.debug(f'starting {option} download')
            with open(outfpath, 'wb') as fp:
                result = ftp.retrbinary(f'RETR {ftp_fpath}', fp.write)
                success = True if result.startswith('226 Transfer complete') else False
            if success:
                config.set(section, option, outfpath)
            logger.debug(f'{option} download successful : {success}')
        ftp.quit()
    config = extract_taxdump(config, dryrun)
    return format_nr(config, dryrun, nproc)

def update_markers(config, options, dryrun):
    """Short summary.

    Parameters
    ----------
    config : type
        Description of parameter `config`.
    options : type
        Description of parameter `options`.
    dryrun : type
        Description of parameter `dryrun`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    section = 'markers'
    if not config.has_section('databases'):
        config.add_section('databases')
    if not config.has_option('databases',section):
        outdir = DEFAULT_CONFIG.get('databases', section)
        config.set('databases', section, outdir)
    else:
        outdir = config.get('databases',section)
    for option in options:
        url = DEFAULT_CONFIG.get('database_urls', option)
        if config.has_option(section, option):
            outfpath = config.get(section, option)
        else:
            outfname = os.path.basename(url)
            outfpath = os.path.join(outdir, outfname)
        if dryrun:
            logger.debug(f'update {section} : {option} : {outfpath}')
            config.set(section, option, outfpath)
            continue
        with requests.Session() as session:
            resp = session.get(url)
        if not resp.ok:
            logger.warning(f'Failed to retrieve {url}')
            continue
        with open(outfpath, 'w') as outfh:
            outfh.write(resp.text)
        config.set(section, option, outfpath)
    return config

def validate_fpaths(config, section):
    """Check all files from section exist and are not empty.

    Parameters
    ----------
    config : type
        Description of parameter `config`.
    section : type
        Description of parameter `section`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    for opt in config.options(section):
        if opt not in DB_SECTIONS.get(section):
            continue
        fp = config.get(section,opt)
        if not os.path.exists(fp) or os.stat(fp).st_size == 0:
            logger.warning(f'removing invalid filepath {fp} : {section} : {opt}')
            config.remove_option(section, opt)
    return config

def update_missing(config, section, dryrun, options=None, nproc=2):
    """Download databases using provided `options` in `section`. If `options` is
    None, all `options` in `section` will be downloaded and formatted.

    Parameters
    ----------
    config: configparser.ConfigParser
        Description of parameter `config`.
    section : str
        Description of parameter `section` (the default is None).
    options : iterable
        Description of parameter `options` (the default is None).
    dryrun : bool
        Description of parameter `dryrun`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    KeyError
        provided `section` is not in DB_SECTIONS

    """
    if section not in DB_SECTIONS:
        raise KeyError(f'section not in DB_SECTIONS : {section}')
    options = set(options) if options else set(DB_SECTIONS.get(section))
    if section == 'ncbi':
        if 'nodes' in options or 'names' in options:
            options.discard('nodes')
            options.discard('names')
            options.add('taxdump')
        config = update_ncbi(config, options, dryrun, nproc)
    if section == 'markers':
        config = update_markers(config, options, dryrun)
    return config

def check_format(config, dryrun, nproc=2):
    """Short summary.

    Parameters
    ----------
    config : type
        Description of parameter `config`.
    dryrun : bool
        Description of parameter `dryrun`.

    Returns
    -------
    type
        Description of returned object.

    Raises
    -------
    ExceptionName
        Why the exception is raised.

    """
    for section,options in DB_SECTIONS.items():
        if not config.has_section(section):
            logger.warning(f'Missing section : {section}')
            config.add_section(section)
            config = update_missing(config, section=section, dryrun=dryrun, nproc=nproc)
            continue
        config = validate_fpaths(config, section)
        missing = set(options) - set(config.options(section))
        if missing:
            logger.warning(f'Missing options : {", ".join(missing)}')
            config = update_missing(config, section=section, options=missing, dryrun=dryrun)

    return config

def configure(config=DEFAULT_CONFIG, dryrun=True, nproc=2):
    """Configures database dependencies necessary to run Autometa.

    Parameters
    ----------
    config : configparser.ConfigParser
        (the default is DEFAULT_CONFIG).
    dryrun : bool
        Log configuration actions but do not perform them. (the default is False).

    Returns
    -------
    configparser.ConfigParser
        config with updated options in respective databases sections

    Raises
    -------
    TypeError
        Provided `config` is not the python built-in ConfigParser type.
    TypeError
        `dryrun` can not be interpretted as a boolean.

    """
    if type(config) is not ConfigParser:
        raise TypeError(f'config is not ConfigParser : {type(config)}')
    if type(dryrun) is not bool:
        raise TypeError(f'dryrun must be True or False. type: {type(dryrun)}')
    return check_format(config, dryrun=dryrun, nproc=nproc)

def main(args):
    config = configure(infpath=args.config, dryrun=args.dryrun)
    if not args.out:
        import sys;sys.exit(0)
    put_config(config, args.out)
    logger.debug(f'{args.out} written.')

if __name__ == '__main__':
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('databases config')
    parser.add_argument('config', help='</path/to/input/database.config>')
    parser.add_argument('--out', help='</path/to/output/database.config>')
    parser.add_argument('--dryrun', help='Log configuration actions but do not perform them.',
        action='store_true', default=False)
    args = parser.parse_args()
    main(args)
