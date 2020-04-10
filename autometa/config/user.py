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

AutometaUser configuration class

"""


import logging
import os

import argparse

from autometa import config
from autometa.config import environ
from autometa.common import utilities

from autometa.common.metagenome import Metagenome
from autometa.config.databases import Databases
from autometa.config.project import Project
from autometa.common.utilities import timeit

logger = logging.getLogger(__name__)


class AutometaUser:
    """AutometaUser Class to handle job submissions.

    Parameters
    ----------
    user_config : str
        </path/to/user/default.config
    dryrun : type
        Description of parameter `dryrun` (the default is True).
    nproc : type
        Description of parameter `nproc` (the default is 2).

    config : config.ConfigParser
        user base database/executable configuration

    Methods
    -------
    configure : NoneType
        Configure user databases and environment.
    dryrun: bool
        If True, will check dependencies and exit.
    nproc: int
        Number of processors to use while downloading/formatting databases.

    """

    def __init__(self, user_config=None, dryrun=True, nproc=2):
        self.dryrun= dryrun
        self.nproc = nproc
        self.config = config.get_config(user_config) if user_config else config.DEFAULT_CONFIG
        if not self.config.has_section('common'):
            self.config.add_section('common')
        self.config.set('common','home_dir', config.AUTOMETA_DIR)

        if self.dryrun:
            self.configure()
            return

    def configure(self):
        """Configure user execution environment and databases.

        Returns
        -------
        NoneType

        """
        # Execution env
        self.config, exe_satisfied = environ.configure(self.config)
        if self.dryrun:
            logger.info(f'Executable dependencies satisfied: {exe_satisfied}')
        # Database env
        dbs = Databases(self.config, dryrun=self.dryrun, nproc=self.nproc)
        self.config = dbs.configure()
        if self.dryrun:
            logger.info(f'Database dependencies satisfied: {dbs.satisfied}')
            return

        if not dbs.satisfied:
            raise LookupError('Database dependencies not satisfied!')
        if not exe_satisfied:
            raise LookupError('Executable dependencies not satisfied!')

    def new_project(self, fpath):
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
        dpath = os.path.dirname(fpath)
        if not os.path.exists(dpath):
            os.makedirs(dpath)
        config.put_config(self.config, fpath)
        return Project(fpath)

    def prepare_binning_args(self, fpath):
        """Prepares metagenome binning run using provided `fpath`.

        Notes
        -----
            This method performs a number of configuration checks to ensure the
            binning run will perform without conflicts.

            1. workspace check: Will construct workspace directory if provided does not exist.
            2. Project check: Will configure a new project if project number is not found in workspace directory.
            3. Metagenome check: Will update if existing with edits or resume if existing without edits. Otherwise will add new metagenome to project.


        Example
        -------
        .. code-block:: python

            #:  Generate namespace - mgargs.files.<file>
            mgargs = prepare_binning_args(mg_config)
            #:  Access file from args - mgargs.files.<file>
            mgargs.files.length_filtered
            #:  Access parameter from args - mgargs.parameters.<parameter>
            mgargs.parameters.length_cutoff

        Parameters
        ----------
        fpath : str
            </path/to/metagenome.config>

        Returns
        -------
        argparse.Namespace
            parameters and files parsed from config.

        Raises
        -------
        ExceptionName
            Why the exception is raised.

        """
        # 1. configure user environment
        self.configure()
        # 2. check workspace exists
        mgargs = config.parse_config(fpath)
        workspace = os.path.realpath(mgargs.parameters.workspace)
        if not os.path.exists(workspace):
            os.makedirs(workspace)
        # 3. check project exists
        proj_name = f'project_{mgargs.parameters.project:03d}'
        project_dirpath = os.path.realpath(os.path.join(workspace,proj_name))
        project_config_fp = os.path.join(project_dirpath, 'project.config')
        if not os.path.exists(project_dirpath) or not os.path.exists(project_config_fp):
            project = self.new_project(project_config_fp)
        else:
            project = Project(project_config_fp)
        # 4. check whether existing or new run with metagenome_num
        metagenome = f'metagenome_{mgargs.parameters.metagenome_num:03d}'
        if metagenome not in project.metagenomes:
            mgargs = project.add(fpath)
            project.save()
            return mgargs
        # If resuming existing metagenome run. Check whether config file has changed.
        old_config_fp = project.metagenomes.get(metagenome)
        old_chksum = utilities.get_checksum(old_config_fp)
        new_chksum = utilities.get_checksum(fpath)
        if old_chksum != new_chksum:
            mgargs = project.update(
                metagenome_num=mgargs.parameters.metagenome_num,
                fpath=fpath)
        project.save()
        return mgargs

    @utilities.timeit
    def run_binning(self, mgargs):
        """Run autometa.

        Parameters
        ----------
        mgargs : argparse.Namespace
            metagenome args

        Returns
        -------
        NoneType

        Raises
        -------
        TODO: Need to enumerate all exceptions raised from within binning pipeline.
        I.e. Demarkate new exception (not yet handled) vs. handled exception.
        Subclassing an AutometaException class may be most appropriate use case here.
        """
        mg = Metagenome(
            assembly=mgargs.files.metagenome,
            outdir=mgargs.parameters.outdir,
            nucl_orfs_fpath=mgargs.files.nucleotide_orfs,
            prot_orfs_fpath=mgargs.files.amino_acid_orfs,
            taxonomy_fpath=mgargs.files.taxonomy,
            fwd_reads=mgargs.files.fwd_reads,
            rev_reads=mgargs.files.rev_reads,
            se_reads=mgargs.files.se_reads,
            taxon_method=mgargs.parameters.taxon_method)
        try:
        # Original (raw) file should not be manipulated so return new object
            mg = mg.length_filter(
                out=mgargs.files.length_filtered,
                cutoff=mgargs.parameters.length_cutoff)
        except FileExistsError as err:
            logger.debug(f'{mgargs.files.length_filtered} already exists. Continuing..')
            mg = Metagenome(
                assembly=mgargs.files.length_filtered,
                outdir=mgargs.parameters.outdir,
                nucl_orfs_fpath=mgargs.files.nucleotide_orfs,
                prot_orfs_fpath=mgargs.files.amino_acid_orfs,
                taxonomy_fpath=mgargs.files.taxonomy,
                fwd_reads=mgargs.files.fwd_reads,
                rev_reads=mgargs.files.rev_reads,
                se_reads=mgargs.files.se_reads,
                taxon_method=mgargs.parameters.taxon_method)
        # I.e. asynchronous execution here (work-queue tasks)
        mg.get_kmers(
            kmer_size=mgargs.parameters.kmer_size,
            normalized=mgargs.files.kmer_normalized,
            out=mgargs.files.kmer_counts,
            multiprocess=mgargs.parameters.kmer_multiprocess,
            nproc=mgargs.parameters.cpus,
            force=mgargs.parameters.force)

        coverages = mg.get_coverages(
            out=mgargs.files.coverages,
            from_spades=mgargs.parameters.cov_from_spades,
            sam=mgargs.files.sam,
            bam=mgargs.files.bam,
            lengths=mgargs.files.lengths,
            bed=mgargs.files.bed)
        # Filter by Kingdom
        kingdoms = mg.get_kingdoms(
            ncbi=mgargs.databases.ncbi,
            usepickle=mgargs.parameters.usepickle,
            blast=mgargs.files.blastp,
            hits=mgargs.files.blastp_hits,
            force=mgargs.parameters.force,
            cpus=mgargs.parameters.cpus)

        if not mgargs.parameters.kingdom in kingdoms:
            raise KeyError(f'{mgargs.parameters.kingdom} not recovered in dataset. Recovered: {", ".join(kingdoms.keys())}')

        mag = kingdoms.get(mgargs.parameters.kingdom)
        bins_df = mag.get_binning(
            method=mgargs.parameters.binning_method,
            kmers=mgargs.files.kmer_counts,
            embedded=mgargs.files.kmer_embedded,
            do_pca=mgargs.parameters.do_pca,
            pca_dims=mgargs.parameters.pca_dims,
            embedding_method=mgargs.parameters.embedding_method,
            coverage=coverages,
            domain=mgargs.parameters.kingdom,
            taxonomy=mgargs.files.taxonomy,
            reverse=mgargs.parameters.reversed,
        )
        binning_cols = ['cluster','completeness','purity']
        bins_df[binning_cols].to_csv(
            mgargs.files.binning,
            sep='\t',
            index=True,
            header=True)

def main(args):
    logger.info(args.user)

if __name__ == '__main__':
    #start_parsing
    import argparse
    import logging as logger
    logger.basicConfig(
        format='%(asctime)s : %(name)s : %(levelname)s : %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p',
        level=logger.DEBUG)
    parser = argparse.ArgumentParser('Concise Functional Description of Script')
    parser.add_argument('user', help='</path/to/user.config>')
    args = parser.parse_args()
    #end_parsing
    main(args)
